import os
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse, StreamingResponse
import io
import xenaPython as xena
import pandas as pd
import numpy as np
from functools import lru_cache

app = FastAPI()

TOIL_HOST = "https://toil.xenahubs.net"
EXPR_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENO_DATASET = "TcgaTargetGTEX_phenotype.txt"

SUBTYPE_MAPPINGS = {
    "Breast Invasive Carcinoma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix", "field": "PAM50Call_RNAseq"},
    "Glioblastoma multiforme": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.GBM.sampleMap/GBM_clinicalMatrix", "field": "GeneExp_Subtype"},
    "Brain Lower Grade Glioma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.LGG.sampleMap/LGG_clinicalMatrix", "field": "IDH_status"},
    "Colon adenocarcinoma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.COAD.sampleMap/COAD_clinicalMatrix", "field": "Subtype_mRNA"}, 
    "Uterine Corpus Endometrial Carcinoma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.UCEC.sampleMap/UCEC_clinicalMatrix", "field": "Subtype_mRNA"},
    "Stomach adenocarcinoma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.STAD.sampleMap/STAD_clinicalMatrix", "field": "Subtype_mRNA"},
    "Thyroid carcinoma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.THCA.sampleMap/THCA_clinicalMatrix", "field": "Subtype_mRNA"},
    "Lung squamous cell carcinoma": {"hub": "https://tcga.xenahubs.net", "dataset": "TCGA.LUSC.sampleMap/LUSC_clinicalMatrix", "field": "expression_subtype"},
}

MUTATION_GENES = ["KRAS", "TP53", "EGFR", "BRAF"]

samples_cache = None
pheno_df_cache = None
mutation_cache = {}

def decode_codes(values, code_map):
    out = []
    for v in values:
        try:
            if v is None or (isinstance(v, float) and np.isnan(v)): out.append(None)
            elif isinstance(v, str) and v.strip().lower() in ("nan", ""): out.append(None)
            else: out.append(code_map.get(int(float(v)), None))
        except: out.append(None)
    return out

def initialize_cache():
    global samples_cache, pheno_df_cache, mutation_cache
    if samples_cache is not None and pheno_df_cache is not None: return
    
    print("Initializing global phenotype cache from Xena... (this might take ~20 seconds)")
    samples_cache = xena.dataset_samples(TOIL_HOST, EXPR_DATASET, None)
    pheno_fields = ["_study", "_sample_type", "primary disease or tissue", "detailed_category"]
    pheno_data = xena.dataset_fetch(TOIL_HOST, PHENO_DATASET, samples_cache, pheno_fields)
    codes_raw = xena.field_codes(TOIL_HOST, PHENO_DATASET, pheno_fields)
    
    code_maps = {}
    for entry in codes_raw:
        name = entry["name"]
        code_maps[name] = ({i: v for i, v in enumerate(entry["code"].split("\t"))} if entry.get("code") else {})
        
    df = pd.DataFrame({
        "sample": samples_cache,
        "study": decode_codes(pheno_data[0], code_maps.get("_study", {})),
        "sample_type": decode_codes(pheno_data[1], code_maps.get("_sample_type", {})),
        "primary_disease_tissue": decode_codes(pheno_data[2], code_maps.get("primary disease or tissue", {})),
        "detailed_category": decode_codes(pheno_data[3], code_maps.get("detailed_category", {})),
    })
    
    # 1. Merge Canonical Subtypes
    df["molecular_subtype"] = None
    for disease, mapping in SUBTYPE_MAPPINGS.items():
        print(f"Fetching subtypes for {disease}...")
        try:
            res = xena.dataset_fetch(mapping["hub"], mapping["dataset"], samples_cache, [mapping["field"]])
            field_codes = xena.field_codes(mapping["hub"], mapping["dataset"], [mapping["field"]])
            if field_codes and field_codes[0].get("code"):
                c_map = {i: v for i, v in enumerate(field_codes[0]["code"].split("\t"))}
                subtype_arr = decode_codes(res[0], c_map)
            else:
                subtype_arr = decode_codes(res[0], {})
            
            mask = df["primary_disease_tissue"] == disease
            for i in range(len(df)):
                if mask[i] and subtype_arr[i] and str(subtype_arr[i]).strip() not in ("nan", "None", ""):
                    df.at[i, "molecular_subtype"] = f"{disease} ({subtype_arr[i]})"
        except Exception as e:
            print(f"Failed to fetch {disease}: {e}")

    # 2. Fetch Mutation Arrays (from PanCanAtlas Hub)
    for m in MUTATION_GENES:
        print(f"Fetching mutation status for {m}...")
        try:
            res = xena.dataset_gene_probe_avg("https://pancanatlas.xenahubs.net", "mc3.v0.2.8.PUBLIC.nonsilentGene.xena", samples_cache, [m])
            if res and res[0].get("scores"):
                mutation_cache[m] = res[0]["scores"][0]
        except Exception as e:
            print(f"Failed to fetch {m} mutations: {e}")
            
    # 3. Fetch MSI Status
    df["msi_status"] = None
    for hub, dataset in [("https://tcga.xenahubs.net", "TCGA.COAD.sampleMap/COAD_clinicalMatrix"),
                         ("https://tcga.xenahubs.net", "TCGA.READ.sampleMap/READ_clinicalMatrix")]:
        try:
            res = xena.dataset_fetch(hub, dataset, samples_cache, ["microsatellite_instability"])
            field_codes = xena.field_codes(hub, dataset, ["microsatellite_instability"])
            if field_codes and field_codes[0].get("code"):
                c_map = {i: v for i, v in enumerate(field_codes[0]["code"].split("\t"))}
                msi_arr = decode_codes(res[0], c_map)
                for i in range(len(df)):
                    if msi_arr[i] and str(msi_arr[i]).strip() not in ("nan", "None", ""):
                        df.at[i, "msi_status"] = str(msi_arr[i]).strip()
        except:
            pass
            
    pheno_df_cache = df
    print("Global phenotype and mutation cache initialized successfully.")


def generate_plot_data(*dfs):
    plot_df = pd.concat([d[["label", "source", "log2_tpm"]] for d in dfs], ignore_index=True)
    medians = plot_df.groupby(["label", "source"])["log2_tpm"].median().sort_values(ascending=False).reset_index()
    
    box_data = []
    for label, source in zip(medians["label"], medians["source"]):
        vals = plot_df[plot_df["label"] == label]["log2_tpm"].dropna().values
        if len(vals) < 3: continue 
        box_data.append({"label": f"{label} (n={len(vals)})", "source": source, "values": vals.tolist()})
        
    summary = (
        plot_df.groupby(["source", "label"])["log2_tpm"]
        .agg(n="count", median="median", mean="mean", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
        .round(3).reset_index().sort_values(["source", "median"], ascending=[True, False])
    )
    return box_data, summary.to_dict(orient="records")


@lru_cache(maxsize=128)
def get_gene_expression(gene: str):
    gene = gene.strip().upper()
    initialize_cache()
    
    print(f"Fetching expression for {gene}...")
    results = xena.dataset_gene_probe_avg(TOIL_HOST, EXPR_DATASET, samples_cache, [gene])
    if not results or not results[0].get("scores") or not results[0]["scores"][0]:
        raise ValueError(f"No expression data returned for gene '{gene}'.")
        
    df = pheno_df_cache.copy()
    df["log2_tpm"] = results[0]["scores"][0]
    
    gtex = df[(df["study"] == "GTEX") & (df["sample_type"] == "Normal Tissue") & (df["detailed_category"].notna())].copy()
    gtex["label"] = gtex["detailed_category"]
    gtex["source"] = "GTEx Normal"

    tcga_base = df[(df["study"] == "TCGA") & (df["sample_type"] == "Primary Tumor") & (df["primary_disease_tissue"].notna())].copy()
    tcga_base["source"] = "TCGA Tumor"
    
    tcga_normal = df[(df["study"] == "TCGA") & (df["sample_type"].astype(str).str.lower().str.contains("normal|control")) & (df["primary_disease_tissue"].notna())].copy()
    tcga_normal["source"] = "TCGA Normal"
    
    # 1. Broad Grouping
    tcga_broad = tcga_base.copy()
    tcga_broad["label"] = tcga_broad["primary_disease_tissue"]
    
    tcga_normal_broad = tcga_normal.copy()
    tcga_normal_broad["label"] = tcga_normal_broad["primary_disease_tissue"] + " (Adj Normal)"
    
    # 2. Canonical Subtypes
    tcga_sub = tcga_base.copy()
    tcga_sub["label"] = tcga_sub["molecular_subtype"].fillna(tcga_sub["primary_disease_tissue"])
    
    tcga_normal_sub = tcga_normal.copy()
    tcga_normal_sub["label"] = tcga_normal_sub["molecular_subtype"].fillna(tcga_normal_sub["primary_disease_tissue"]) + " (Adj Normal)"
    
    broad_plot, broad_summary = generate_plot_data(tcga_broad, gtex, tcga_normal_broad)
    sub_plot, sub_summary = generate_plot_data(tcga_sub, gtex, tcga_normal_sub)
    
    payload = {"gene": gene, "broad": {"plot_data": broad_plot, "summary": broad_summary}, "subtype": {"plot_data": sub_plot, "summary": sub_summary}}
    
    # 3. Mutational Subtypes
    for m in MUTATION_GENES:
        if m in mutation_cache and mutation_cache[m]:
            mut_scores = mutation_cache[m]
            
            # Apply to TCGA base cohort
            df_mut = tcga_base.copy()
            df_mut["mut_stat"] = [mut_scores[idx] for idx in df_mut.index]
            df_mut = df_mut[df_mut["mut_stat"].notna()]
            
            labels = []
            for stat, disease in zip(df_mut["mut_stat"], df_mut["primary_disease_tissue"]):
                try:
                    is_mut = float(stat) > 0.0
                except (ValueError, TypeError):
                    is_mut = False
                    
                if is_mut: labels.append(f"{disease} ({m} Mutant)")
                else: labels.append(f"{disease} ({m} WT)")
                
            df_mut["label"] = labels
            
            m_plot, m_sum = generate_plot_data(df_mut, gtex, tcga_normal_broad)
            if m_plot:
                payload[f"mut_{m}"] = {"plot_data": m_plot, "summary": m_sum}
                
    # 4. MSI Status
    if tcga_base["msi_status"].notna().any():
        df_msi = tcga_base.copy()
        labels = []
        for disease, msi in zip(df_msi["primary_disease_tissue"], df_msi["msi_status"]):
            if pd.notna(msi) and str(msi).strip() not in ("nan", "None", ""):
                labels.append(f"{disease} (MSI: {msi})")
            else:
                labels.append(disease)
        df_msi["label"] = labels
        
        msi_plot, msi_sum = generate_plot_data(df_msi, gtex, tcga_normal_broad)
        if msi_plot:
            payload["msi_status"] = {"plot_data": msi_plot, "summary": msi_sum}
            
    return payload

@app.get("/api/expression")
def api_expression(gene: str):
    try:
        return JSONResponse(content=get_gene_expression(gene))
    except ValueError as e: return JSONResponse(status_code=404, content={"detail": str(e)})
    except Exception as e: return JSONResponse(status_code=500, content={"detail": f"Internal error: {e}"})

@app.get("/api/top_targets")
def api_top_targets():
    csv_path = "surface_selectivity_scores.csv"
    if not os.path.exists(csv_path):
        return JSONResponse(status_code=404, content={"detail": "Selectivity scores not found."})
    df = pd.read_csv(csv_path)
    df = df.replace({np.nan: None})
    return JSONResponse(content=df.to_dict(orient="records"))

@app.get("/api/download_csv")
def api_download_csv():
    csv_path = "surface_selectivity_scores.csv"
    if not os.path.exists(csv_path):
        return JSONResponse(status_code=404, content={"detail": "CSV not found."})
        
    df = pd.read_csv(csv_path)
    if os.path.exists("surface_genes.txt"):
        with open("surface_genes.txt", "r") as f:
            surface_genes = set([line.strip() for line in f if line.strip()])
        df["Is_Surface_Protein"] = df["Gene"].apply(lambda x: "Yes" if x in surface_genes else "No")
    
    stream = io.StringIO()
    df.to_csv(stream, index=False)
    response = StreamingResponse(iter([stream.getvalue()]), media_type="text/csv")
    response.headers["Content-Disposition"] = "attachment; filename=surface_selectivity_scores.csv"
    return response

@app.get("/api/surface_genes")
def api_surface_genes():
    if not os.path.exists("surface_genes.txt"):
        return JSONResponse(status_code=404, content={"detail": "surface_genes.txt not found."})
    with open("surface_genes.txt", "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    return JSONResponse(content={"genes": genes})

static_dir = os.path.join(os.path.dirname(__file__), "static")
os.makedirs(static_dir, exist_ok=True)
app.mount("/", StaticFiles(directory=static_dir, html=True), name="static")

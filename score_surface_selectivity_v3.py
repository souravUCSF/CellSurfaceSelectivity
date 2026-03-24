import os
import gzip
import urllib.request
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PHENOTYPE_DATAFILE = "data/phenotype.txt.gz"
EXPR_DATAFILE = "data/TcgaTargetGtex_rsem_gene_tpm.gz"
PROBEMAP_FILE = "data/probeMap.tsv"

def fetch_probemap():
    mapper = {}
    with open(PROBEMAP_FILE, "r") as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                mapper[parts[0].split('.')[0]] = parts[1]
    return mapper

def get_phenotype_annotations():
    if not os.path.exists(PHENOTYPE_DATAFILE):
        urllib.request.urlretrieve("https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz", PHENOTYPE_DATAFILE)
        
    df = pd.read_csv(PHENOTYPE_DATAFILE, sep='\t', compression='gzip', encoding='latin1')
    df.rename(columns={"primary disease or tissue": "primary_disease"}, inplace=True)
    df["group"] = "Unknown"
    
    for i, row in df.iterrows():
        study = str(row["_study"]).strip().upper()
        sample_type = str(row["_sample_type"]).strip().lower()
        disease = str(row["primary_disease"]).strip()
        cat = str(row["detailed_category"]).strip()
        
        if disease == "" or disease.lower() == "nan":
            continue
            
        # Match exactly what app.py boxplots use
        if study == "TCGA" and sample_type == "primary tumor":
            df.at[i, "group"] = f"Tumor: {disease}"
        elif study == "TCGA" and ("normal" in sample_type or "control" in sample_type) and disease != "" and disease.lower() != "nan":
            df.at[i, "group"] = f"TCGANormal: {disease}"
        elif study == "GTEX" and sample_type == "normal tissue" and cat != "" and cat.lower() != "nan":
            df.at[i, "group"] = f"Normal: {cat}"
            
    return df[df["group"] != "Unknown"]

def calculate_tau(expression_profile):
    if len(expression_profile) <= 1: return 0.0
    # Data is log2(TPM + 0.001), convert back to linear TPM + 0.001 to assure positivity for Tau
    tpm_vals = (2 ** expression_profile)
    max_expr = np.max(tpm_vals)
    if max_expr == 0: return 0.0
    return np.sum(1 - (tpm_vals / max_expr)) / (len(tpm_vals) - 1)

def run_transcriptome_pandas(mapper, pheno_df):
    sample_to_group = dict(zip(pheno_df["sample"], pheno_df["group"]))
    
    with gzip.open(EXPR_DATAFILE, "rt") as f:
        header_samples = f.readline().strip().split('\t')[1:]
    
    col_to_group = []
    usecols = [0]
    
    for idx, sample in enumerate(header_samples):
        if sample in sample_to_group:
            col_to_group.append(sample_to_group[sample])
            usecols.append(idx + 1)
            
    tumor_groups = sorted(list(set(g for g in col_to_group if g.startswith("Tumor:"))))
    normal_groups = sorted(list(set(g for g in col_to_group if g.startswith("Normal:"))))
    tcga_normal_groups = sorted(list(set(g for g in col_to_group if g.startswith("TCGANormal:"))))
    
    print(f"Mapped {len(usecols)-1} samples to {len(tumor_groups)} tumor, {len(normal_groups)} GTEx, and {len(tcga_normal_groups)} TCGA normal groups.")
    if len(tumor_groups) == 0:
        return pd.DataFrame()
        
    chunksize = 200 # Minimal chunksize to avoid memory swapping
    results = []
    
    chunk_iter = pd.read_csv(EXPR_DATAFILE, sep='\t', compression='gzip', usecols=usecols, chunksize=chunksize, index_col=0)
    
    processed = 0
    for chunk in chunk_iter:
        chunk.columns = col_to_group
        medians_df = chunk.T.groupby(chunk.columns).median().T
        
        tumor_medians = medians_df[tumor_groups]
        normal_medians = medians_df[normal_groups]
        tcga_n_medians = medians_df[tcga_normal_groups]
        normal_mask = chunk.columns.str.startswith("Normal:")
        
        for gene_id, row in chunk.iterrows():
            gene_id_str = str(gene_id).split(".")[0]
            symbol = mapper.get(gene_id_str, gene_id_str)
            
            t_meds = tumor_medians.loc[gene_id].values
            n_meds = normal_medians.loc[gene_id].values
            tcga_n_meds = tcga_n_medians.loc[gene_id].values
            
            t_max_idx = np.argmax(t_meds)
            t_max = t_meds[t_max_idx]
            max_tumor_group = tumor_groups[t_max_idx]
            
            n_max = np.max(n_meds) if len(n_meds) > 0 else 0.0
            t_adjn_max = np.max(tcga_n_meds) if len(tcga_n_meds) > 0 else 0.0
            
            delta_tn = max(0, t_max - n_max)
            delta_t_adjn = max(0, t_max - t_adjn_max)
            
            tau = calculate_tau(t_meds)
            
            score = delta_tn * tau
            score_t_adjn = delta_t_adjn * tau
            
            pval = 1.0
            if delta_tn > 0 and tau > 0:
                top_tumor_raw = row[chunk.columns == max_tumor_group].values
                all_normal_vals = row[normal_mask].values
                
                top_tumor_raw = top_tumor_raw[pd.notna(top_tumor_raw)]
                all_normal_vals = all_normal_vals[pd.notna(all_normal_vals)]
                
                if len(top_tumor_raw) > 0 and len(all_normal_vals) > 0:
                    try:
                        _, pval = mannwhitneyu(top_tumor_raw, all_normal_vals, alternative='greater')
                    except ValueError:
                        pass
                        
            results.append({
                "Gene": symbol, "Score": score, "p_value": pval,
                "Tumor_Tau": tau, "Delta_TN": delta_tn, "T_max": t_max,
                "max_tumor_type": max_tumor_group.replace("Tumor: ", ""),
                "N_max": n_max, "ensembl_id": gene_id_str,
                "Score_T_AdjN": score_t_adjn, "TAdjN_max": t_adjn_max
            })
            
        processed += len(chunk)
        if processed % 1000 == 0:
            print(f"Processed {processed} genes.", flush=True)

    df = pd.DataFrame(results)
    
    if len(df) > 0:
        from scipy.stats import false_discovery_control
        df["p_adj"] = false_discovery_control(df["p_value"])
        cols = ["Gene", "Score", "Score_T_AdjN", "p_value", "p_adj", "Tumor_Tau", "Delta_TN", "T_max", "max_tumor_type", "N_max", "TAdjN_max", "ensembl_id"]
        df = df[cols].sort_values("Score", ascending=False)
        
    return df

def main():
    mapper = fetch_probemap()
    pheno_df = get_phenotype_annotations()
    final_scores = run_transcriptome_pandas(mapper, pheno_df)
    
    if len(final_scores) > 0:
        final_scores.to_csv("surface_selectivity_scores.csv", index=False)
        print("\nTop 15 Selectivity Genes:")
        print(final_scores.head(15)[["Gene", "Score", "p_adj", "max_tumor_type", "T_max", "N_max"]].to_string(index=False))
    else:
        print("Error: No scores were generated!")

if __name__ == "__main__":
    main()

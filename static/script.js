const form = document.getElementById('search-form');
const input = document.getElementById('gene-input');
const resultsContainer = document.getElementById('results-container');
const spinner = document.getElementById('loading-spinner');
const errorMsg = document.getElementById('error-message');
const plotContainer = document.getElementById('plot-container');
const plotTitle = document.getElementById('plot-title');
const toggleContainer = document.getElementById('view-toggle-container');
const groupSelect = document.getElementById('group-select');
const topStatsTable = document.getElementById('top-stats-table').querySelector('tbody');
const bottomStatsTable = document.getElementById('bottom-stats-table').querySelector('tbody');

const TCGA_COLOR = 'rgba(239, 68, 68, 0.75)'; // #ef4444
const GTEX_COLOR = 'rgba(59, 130, 246, 0.75)'; // #3b82f6
const LINE_COLOR = 'rgba(255, 255, 255, 0.8)';

// Selectivity Table Elements
const targetsContainer = document.getElementById('targets-container');
const mainTargetsTable = document.getElementById('main-targets-table').querySelector('tbody');
const filterGene = document.getElementById('filter-gene');
const filterTumor = document.getElementById('filter-tumor');
const filterSurface = document.getElementById('filter-surface-only');

const filterScore = document.getElementById('filter-score');
const filterTau = document.getElementById('filter-tau');
const filterPval = document.getElementById('filter-pval');
const filterTmax = document.getElementById('filter-tmax');
const filterNmax = document.getElementById('filter-nmax');
const filterDelta = document.getElementById('filter-delta');

const valScore = document.getElementById('val-score');
const valTau = document.getElementById('val-tau');
const valPval = document.getElementById('val-pval');
const valTmax = document.getElementById('val-tmax');
const valNmax = document.getElementById('val-nmax');
const valDelta = document.getElementById('val-delta');

const backToTableBtn = document.getElementById('back-to-table-btn');
const pagePrevBtn = document.getElementById('page-prev');
const pageNextBtn = document.getElementById('page-next');
const pageInfo = document.getElementById('page-info');
const paginationControls = document.getElementById('pagination-controls');

let currentDataPayload = null;
let currentGene = "";

// Table tracking state
let allTargets = [];
let filteredTargets = [];
let surfaceGenes = new Set();
let currentPage = 1;
const rowsPerPage = 50;

let sortColumn = 'Score';
let sortDirection = -1; // -1 for desc, 1 for asc

function sortData() {
    filteredTargets.sort((a, b) => {
        let valA = a[sortColumn];
        let valB = b[sortColumn];
        
        if (valA !== null && valA !== undefined && valB !== null && valB !== undefined && !isNaN(valA) && !isNaN(valB)) {
            valA = Number(valA);
            valB = Number(valB);
        } else {
            valA = String(valA || "").toLowerCase();
            valB = String(valB || "").toLowerCase();
        }
        
        if (valA < valB) return -1 * sortDirection;
        if (valA > valB) return 1 * sortDirection;
        return 0;
    });
}

// Initialize on load
document.addEventListener('DOMContentLoaded', async () => {
    try {
        const [resTargets, resSurface] = await Promise.all([
            fetch('/api/top_targets'),
            fetch('/api/surface_genes')
        ]);
        
        if (resSurface.ok) {
            const sd = await resSurface.json();
            sd.genes.forEach(g => surfaceGenes.add(g.toUpperCase()));
        }

        if (resTargets.ok) {
            allTargets = await resTargets.json();
            const defaultTh = document.querySelector('th[data-sort="Score"]');
            if (defaultTh) defaultTh.classList.add('sort-desc');
            
            applyFilters();
        }
    } catch(e) {
        console.error("Failed to load initial data", e);
    }
});

function renderTablePage() {
    mainTargetsTable.innerHTML = '';
    const totalPages = Math.ceil(filteredTargets.length / rowsPerPage) || 1;
    
    if (currentPage > totalPages) currentPage = totalPages;
    if (currentPage < 1) currentPage = 1;
    
    const startIdx = (currentPage - 1) * rowsPerPage;
    const endIdx = startIdx + rowsPerPage;
    const rowsToShow = filteredTargets.slice(startIdx, endIdx);
    
    rowsToShow.forEach(target => {
        const tr = document.createElement('tr');
        // Format p-value
        let pvalStr = Number(target.p_adj).toExponential(2);
        if (Number(target.p_adj) < 1e-300) pvalStr = "< 1e-300";
        else if (Number(target.p_adj) > 0.01) pvalStr = Number(target.p_adj).toFixed(3);
        
        tr.innerHTML = `
            <td style="font-weight: 600; color: #60a5fa">${target.Gene}</td>
            <td style="font-weight: 600">${Number(target.Score).toFixed(2)}</td>
            <td style="font-weight: 600; color: #a8b8d8">${Number(target.Score_T_AdjN || 0).toFixed(2)}</td>
            <td>${target.max_tumor_type}</td>
            <td style="color: #cbd5e1">${pvalStr}</td>
            <td>${Number(target.Tumor_Tau).toFixed(2)}</td>
            <td>${Number(target.T_max).toFixed(2)}</td>
            <td>${Number(target.N_max).toFixed(2)}</td>
            <td>${Number(target.Delta_TN).toFixed(2)}</td>
        `;
        tr.addEventListener('click', () => {
            input.value = target.Gene;
            form.dispatchEvent(new Event('submit'));
        });
        mainTargetsTable.appendChild(tr);
    });
    
    pageInfo.textContent = `Page ${currentPage} of ${totalPages} (${filteredTargets.length} results)`;
    pagePrevBtn.disabled = currentPage === 1;
    pageNextBtn.disabled = currentPage === totalPages;
    
    if (totalPages > 1) {
        paginationControls.classList.remove('hidden');
    } else {
        paginationControls.classList.add('hidden');
    }
}

function applyFilters() {
    const geneQuery = filterGene.value.trim().toUpperCase();
    const tumorQuery = filterTumor.value.trim().toLowerCase();
    const surfaceOnly = filterSurface.checked;
    
    // Numeric values
    const minScore = parseFloat(filterScore.value);
    const minTau = parseFloat(filterTau.value) || 0;
    const maxPval = parseFloat(filterPval.value) || 1;
    const minTmax = parseFloat(filterTmax.value) || -10;
    const maxNmax = parseFloat(filterNmax.value) || 15;
    const minDelta = parseFloat(filterDelta.value) || 0;
    
    filteredTargets = allTargets.filter(t => {
        const matchGene = !geneQuery || String(t.Gene).toUpperCase().includes(geneQuery);
        const matchTumor = !tumorQuery || String(t.max_tumor_type).toLowerCase().includes(tumorQuery);
        const matchSurface = !surfaceOnly || surfaceGenes.has(String(t.Gene).toUpperCase()) || surfaceGenes.has(t.ensembl_id);
        
        const matchNum = (t.Score >= minScore) &&
                         (t.Tumor_Tau >= minTau) &&
                         (t.p_adj <= maxPval) &&
                         (t.T_max >= minTmax) &&
                         (t.N_max <= maxNmax) &&
                         (t.Delta_TN >= minDelta);
                         
        return matchGene && matchTumor && matchSurface && matchNum;
    });
    
    sortData();
    currentPage = 1;
    renderTablePage();
    renderScatterPlot();
}

document.querySelectorAll('th.sortable').forEach(th => {
    th.addEventListener('click', () => {
        const column = th.dataset.sort;
        if (sortColumn === column) {
            sortDirection *= -1;
        } else {
            sortColumn = column;
            sortDirection = -1;
        }
        
        document.querySelectorAll('th.sortable').forEach(header => {
            header.classList.remove('sort-asc', 'sort-desc');
        });
        
        if (sortDirection === 1) {
            th.classList.add('sort-asc');
        } else {
            th.classList.add('sort-desc');
        }
        
        sortData();
        currentPage = 1;
        renderTablePage();
    });
});

filterGene.addEventListener('input', applyFilters);
filterTumor.addEventListener('input', applyFilters);
filterSurface.addEventListener('change', applyFilters);

[filterScore, filterTau, filterPval, filterTmax, filterNmax, filterDelta].forEach(input => {
    if (!input) return;
    input.addEventListener('input', (e) => {
        const targetLabel = document.getElementById('val-' + e.target.id.split('-')[1]);
        if(targetLabel) targetLabel.textContent = e.target.value;
        applyFilters();
    });
});

pagePrevBtn.addEventListener('click', () => {
    if (currentPage > 1) {
        currentPage--;
        renderTablePage();
    }
});

const downloadCsvBtn = document.getElementById('download-csv-btn');
if (downloadCsvBtn) {
    downloadCsvBtn.addEventListener('click', (e) => {
        e.preventDefault();
        window.location.href = `/api/download_csv?t=${new Date().getTime()}`;
    });
}

function renderScatterPlot() {
    const plotDiv = document.getElementById('targets-plotly-div');
    if (!plotDiv) return;

    const xData = [];
    const yData = [];
    const textData = [];

    filteredTargets.forEach(target => {
        if (target.Tumor_Tau !== null && target.Tumor_Tau !== undefined && 
            target.Delta_TN !== null && target.Delta_TN !== undefined) {
            
            xData.push(Number(target.Tumor_Tau));
            yData.push(Number(target.Delta_TN));
            textData.push(`${target.Gene}<br>Max Tumor: ${target.max_tumor_type}`);
        }
    });

    const trace = {
        x: xData,
        y: yData,
        text: textData,
        mode: 'markers',
        type: 'scattergl',
        marker: {
            size: 6,
            color: 'rgba(99, 102, 241, 0.6)',
            line: { color: 'rgba(255, 255, 255, 0.8)', width: 0.5 }
        },
        hovertemplate: '<b>%{text}</b><br>Tau: %{x:.2f}<br>Delta_TN: %{y:.2f}<extra></extra>'
    };

    const layout = {
        paper_bgcolor: 'transparent',
        plot_bgcolor: 'transparent',
        margin: { l: 60, r: 30, t: 20, b: 50 },
        xaxis: {
            title: 'Tumor Tau (Specificity)',
            gridcolor: 'rgba(255, 255, 255, 0.08)',
            zerolinecolor: 'rgba(255, 255, 255, 0.2)',
            tickfont: { color: '#94a3b8', family: 'Inter' },
            titlefont: { color: '#e2e8f0', family: 'Inter', size: 13 }
        },
        yaxis: {
            title: 'Delta_TN (Expression Difference)',
            gridcolor: 'rgba(255, 255, 255, 0.08)',
            zerolinecolor: 'rgba(255, 255, 255, 0.2)',
            tickfont: { color: '#94a3b8', family: 'Inter' },
            titlefont: { color: '#e2e8f0', family: 'Inter', size: 13 }
        },
        hovermode: 'closest',
        font: { family: 'Inter, sans-serif' }
    };

    Plotly.react(plotDiv, [trace], layout, { responsive: true, displaylogo: false });
}

pageNextBtn.addEventListener('click', () => {
    const totalPages = Math.ceil(filteredTargets.length / rowsPerPage);
    if (currentPage < totalPages) {
        currentPage++;
        renderTablePage();
    }
});

backToTableBtn.addEventListener('click', () => {
    resultsContainer.classList.add('hidden');
    targetsContainer.classList.remove('hidden');
    targetsContainer.scrollIntoView({ behavior: 'smooth', block: 'start' });
});

form.addEventListener('submit', async (e) => {
    e.preventDefault();
    const gene = input.value.trim().toUpperCase();
    if (!gene) return;
    
    currentGene = gene;
    
    targetsContainer.classList.add('hidden');
    resultsContainer.classList.remove('hidden');
    plotContainer.classList.add('hidden');
    errorMsg.classList.add('hidden');
    spinner.classList.remove('hidden');
    toggleContainer.classList.add('hidden');
    
    resultsContainer.scrollIntoView({ behavior: 'smooth', block: 'start' });
    
    try {
        const response = await fetch(`/api/expression?gene=${encodeURIComponent(gene)}`);
        const data = await response.json();
        
        if (!response.ok) {
            throw new Error(data.detail || 'Failed to fetch expression data');
        }
        
        currentDataPayload = data;
        
        // Reset dropdown to standard
        groupSelect.value = "broad";
        
        toggleContainer.classList.remove('hidden');
        renderActiveView();
        
    } catch (err) {
        spinner.classList.add('hidden');
        errorMsg.textContent = err.message;
        errorMsg.classList.remove('hidden');
    }
});

groupSelect.addEventListener('change', () => {
    if (currentDataPayload) {
        renderActiveView();
    }
});

function renderActiveView() {
    const key = groupSelect.value;
    const viewData = currentDataPayload[key];
    
    if (!viewData) {
        // Safe fallback if mutation datasets were down during API startup
        plotContainer.classList.add('hidden');
        errorMsg.textContent = "Detailed subgroup mappings are unavailable for this selection dynamically.";
        errorMsg.classList.remove('hidden');
        return;
    }
    
    errorMsg.classList.add('hidden');
    spinner.classList.add('hidden');
    plotContainer.classList.remove('hidden');
    
    plotTitle.textContent = `${currentGene} Expression Profile`;
    
    // ----- Render Plotly Chart -----
    const traces = [];
    viewData.plot_data.forEach(item => {
        traces.push({
            type: 'box',
            x: item.values,
            name: item.label,
            orientation: 'h',
            boxpoints: false,
            line: { color: LINE_COLOR, width: 1.5 },
            fillcolor: item.source === 'TCGA Tumor' ? TCGA_COLOR : (item.source === 'TCGA Normal' ? 'rgba(16, 185, 129, 0.75)' : GTEX_COLOR),
            marker: { color: 'transparent' },
            hovertemplate: `<b>%{name}</b><br>Median: %{x:.2f}<br>Min: %{min:.2f} · Max: %{max:.2f}<extra></extra>`
        });
    });
    
    const plotHeight = Math.max(700, viewData.plot_data.length * 22);
    
    const layout = {
        title: false,
        paper_bgcolor: 'transparent',
        plot_bgcolor: 'transparent',
        margin: { l: 280, r: 30, t: 20, b: 60 },
        height: plotHeight,
        xaxis: {
            title: 'log₂(TPM + 0.001)',
            gridcolor: 'rgba(255, 255, 255, 0.08)',
            zerolinecolor: 'rgba(255, 255, 255, 0.2)',
            tickfont: { color: '#94a3b8', family: 'Inter' },
            titlefont: { color: '#e2e8f0', family: 'Inter', size: 14 }
        },
        yaxis: {
            gridcolor: 'rgba(255, 255, 255, 0)',
            tickfont: { color: '#e2e8f0', family: 'Inter', size: 11 },
            autorange: 'reversed'
        },
        showlegend: false,
        hovermode: 'closest',
        font: { family: 'Inter, sans-serif' }
    };
    
    const config = { responsive: true, displayModeBar: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d'] };
    Plotly.react('plotly-div', traces, layout, config);
    
    // ----- Render Summary Tables -----
    const globalSorted = [...viewData.summary].sort((a, b) => b.median - a.median);
    const top10 = globalSorted.slice(0, 10);
    const bottom10 = globalSorted.slice(-10).reverse();
    
    function populateTable(tbody, dataArr) {
        tbody.innerHTML = '';
        dataArr.forEach(row => {
            const tr = document.createElement('tr');
            const badgeColor = row.source === 'TCGA Tumor' ? 'rgba(239, 68, 68, 0.2)' : (row.source === 'TCGA Normal' ? 'rgba(16, 185, 129, 0.2)' : 'rgba(59, 130, 246, 0.2)');
            const textColor = row.source === 'TCGA Tumor' ? '#fca5a5' : (row.source === 'TCGA Normal' ? '#6ee7b7' : '#93c5fd');
            
            tr.innerHTML = `
                <td style="font-weight: 500">${row.label}</td>
                <td><span class="badge" style="background: ${badgeColor}; color: ${textColor}; border: 1px solid ${textColor}">${row.source}</span></td>
                <td>${row.n.toLocaleString()}</td>
                <td style="font-weight: 600; color: #fff">${row.median.toFixed(2)}</td>
            `;
            tbody.appendChild(tr);
        });
    }
    
    populateTable(topStatsTable, top10);
    populateTable(bottomStatsTable, bottom10);
}

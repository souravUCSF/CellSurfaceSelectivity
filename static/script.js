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

// Basic UI tracking state
let currentDataPayload = null;
let currentGene = "";

document.addEventListener('DOMContentLoaded', () => {
    // Intentionally simplified public deployment mode
});

form.addEventListener('submit', async (e) => {
    e.preventDefault();
    const gene = input.value.trim().toUpperCase();
    if (!gene) return;
    
    currentGene = gene;
    
    // Switch to results view automatically since table is gone
    resultsContainer.classList.remove('hidden');
    
    plotContainer.classList.add('hidden');
    errorMsg.classList.add('hidden');
    spinner.classList.remove('hidden');
    toggleContainer.classList.add('hidden');
    
    try {
        const response = await fetch(`/api/expression?gene=${gene}`);
        if (!response.ok) {
            const err = await response.json();
            throw new Error(err.detail || 'Failed to fetch data');
        }
        
        currentDataPayload = await response.json();
        renderActiveView();
        
    } catch (err) {
        errorMsg.textContent = err.message;
        errorMsg.classList.remove('hidden');
    } finally {
        spinner.classList.add('hidden');
    }
});

groupSelect.addEventListener('change', renderActiveView);

function renderActiveView() {
    if (!currentDataPayload) return;
    
    const selectedGroup = groupSelect.value;
    const viewData = currentDataPayload[selectedGroup];
    
    if (!viewData || !viewData.plot_data || viewData.plot_data.length === 0) {
        // Fallback or hide
        document.getElementById('plotly-div').innerHTML = `<p style="text-align:center; padding: 2rem; color: #94a3b8;">Data not available for this subset.</p>`;
        topStatsTable.innerHTML = '';
        bottomStatsTable.innerHTML = '';
        plotContainer.classList.remove('hidden');
        toggleContainer.classList.remove('hidden');
        plotTitle.textContent = `${currentGene} Expression`;
        return;
    }
    
    const dataArr = viewData.plot_data;
    
    // Prepare Plotly Box Traces
    const traces = dataArr.map(item => {
        return {
            x: item.values,
            type: 'box',
            name: item.label,
            orientation: 'h',
            boxpoints: false,
            line: { color: LINE_COLOR, width: 1.5 },
            fillcolor: item.source === 'TCGA Tumor' ? TCGA_COLOR : (item.source === 'TCGA Normal' ? 'rgba(16, 185, 129, 0.75)' : GTEX_COLOR),
            marker: { color: 'transparent' },
            hovertemplate: `<b>%{name}</b><br>Median: %{x:.2f}<br>Min: %{min:.2f} · Max: %{max:.2f}<extra></extra>`
        };
    });

    const layout = {
        height: Math.max(800, dataArr.length * 28 + 120),
        margin: { l: 300, r: 20, t: 40, b: 60 },
        xaxis: { title: 'Log2(TPM + 1)', color: '#e2e8f0', gridcolor: 'rgba(255,255,255,0.05)', linecolor: 'rgba(255,255,255,0.1)' },
        yaxis: { tickfont: { size: 12, color: '#f8fafc' }, fixedrange: true, autorange: 'reversed' },
        paper_bgcolor: 'transparent',
        plot_bgcolor: 'transparent',
        showlegend: false,
        font: { family: 'Inter, sans-serif' }
    };

    const config = { responsive: true, displayModeBar: false };
    Plotly.newPlot('plotly-div', traces, layout, config);
    
    // Render Stats Highlights
    const globalSorted = [...viewData.summary].sort((a, b) => b.median - a.median);
    const top10 = globalSorted.slice(0, 10);
    const bottom10 = globalSorted.slice(-10).reverse();
    
    populateStatsTable(topStatsTable, top10);
    populateStatsTable(bottomStatsTable, bottom10);
    
    plotTitle.textContent = `${currentGene} Expression Landscape`;
    plotContainer.classList.remove('hidden');
    toggleContainer.classList.remove('hidden');
}

function populateStatsTable(tbody, dataArr) {
    tbody.innerHTML = '';
    dataArr.forEach(row => {
        const tr = document.createElement('tr');
        const badgeColor = row.source === 'TCGA Tumor' ? 'rgba(239, 68, 68, 0.2)' : (row.source === 'TCGA Normal' ? 'rgba(16, 185, 129, 0.2)' : 'rgba(59, 130, 246, 0.2)');
        const textColor = row.source === 'TCGA Tumor' ? '#fca5a5' : (row.source === 'TCGA Normal' ? '#6ee7b7' : '#93c5fd');
        
        tr.innerHTML = `
            <td style="font-weight: 500">${row.label}</td>
            <td><span class="source-badge" style="background:${badgeColor};color:${textColor}">${row.source}</span></td>
            <td>${row.n.toLocaleString()}</td>
            <td style="font-weight: 600">${Number(row.median).toFixed(2)}</td>
        `;
        tbody.appendChild(tr);
    });
}

const search = document.getElementById('search');

const matchList = document.getElementById('match-list');

// Serach genes.json and filter it

const searchGenes = async searchText => {
    const res = await fetch('static/genes.json');
    const genes = await res.json();

    let matches = genes.filter(gene => {
        const regex = new RegExp(`^${searchText}`, 'gi');
        return regex.test(gene);
    });

    if (searchText.length === 0) {
        matches = [];
        matchList.innerHTML = '';

    }

    outputHtml(matches);
}

const outputHtml = matches => {
    if (matches.length > 0) {
        
        const html = matches.map(match => `
            <div class="card mb-1"> 
                ${match}
            </div>
        `).join('');

    matchList.innerHTML = html;
    }
}

search.addEventListener('input', () => searchGenes(search.value));



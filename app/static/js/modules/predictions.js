export async function get_prediction_studies() {
    var response = await $.ajax({
        url: "api/get_prediction",
        contentType: 'application/json',
        type: "GET",
    })

    const curr_pred = response['curr_prediction'];
    const gse = curr_pred[1].split('-')[0]
    document.getElementById('study1').innerHTML = `<p class="tooltip ml-2" style="margin-block-start: 0em !important; display: inline; font-size: 1rem; left: -4px;">
        <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gse}" target="_blank" rel="noopener noreferrer">${gse}</a>
        <span class="tooltiptext mt-4 text-left" style="position: absolute;">
            Title: ${response['gse_title']}<br><br>
            Signature: ${response['gse_sig']}
        </span></p>`
    document.getElementById('study2').innerHTML = `<p class="tooltip ml-2" style="margin-block-start: 0em !important; display: inline; font-size: 1rem; left: -4px;">
    <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/${curr_pred[4]}" target="_blank"  rel="noopener noreferrer">${curr_pred[4]}</a>
    <span class="tooltiptext mt-4 text-left">
    ${curr_pred[5]}
    </span></p>`
    
    
    

}

export async function get_curr_prediction() {

    const selecter_div = document.getElementById('gpt-hypothesis-selector')
    const hypothesis_div = document.getElementById('gpt-hypothesis')

    selecter_div.innerHTML = 'Loading...'
    var response = await $.ajax({
        url: "api/get_prediction",
        contentType: 'application/json',
        type: "GET"
    })

    const date_options = response['date_options']
    const dates = Object.keys(date_options)
    dates.reverse()
    var selecter = `<div class='text-center'><p class='mb-1 text-bold'>View previous hypotheses: </p><select class="m-2 libpicker" data-style="btn-primary" onchange="get_prediction_date(this)" data-width="500px">`
    for (var i = 0; i < dates.length; i++) {
        var row_num = date_options[dates[i]];
        selecter += `<option value=${row_num}>${dates[i]}</option>`;
    }
    selecter += `</select></div>`

    const curr_pred = response['curr_prediction'];

    const gse = curr_pred[1].split('-')[0]

    const currentUrl = window.location.href
    var gse_gene_viewer = currentUrl.split('/')
    gse_gene_viewer = gse_gene_viewer.filter((x) => x != 'hypotheses')
    gse_gene_viewer.push(gse)
    gse_gene_viewer = gse_gene_viewer.join('/')

    document.getElementById("gpt-hypothesis-title").innerText = `Today's Hypothesis ${dates[0]}`

    const gse_link = `<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gse}" target="_blank"  rel="noopener noreferrer">${gse}</a> `
    const pmc_link = `<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/${curr_pred[4]}" target="_blank"  rel="noopener noreferrer">${curr_pred[4]}</a> `
    var curr_pred_str = "<div class='text-center text-bold'>"
    curr_pred_str += `${gse_link + curr_pred[1].split('-').slice(1).join(' ').replaceAll('_', ' ')} (${curr_pred[2]}),<br> ${pmc_link + curr_pred[3].split('-').slice(1).join(' ').replaceAll('_', ' ')} (${curr_pred[5]}) (${curr_pred[6]})<br>`
    curr_pred_str += `P-value: ${parseFloat(curr_pred[7]).toExponential(2)}, Adj. P-value: ${parseFloat(curr_pred[8]).toExponential(2)}, Odds ratio: ${parseFloat(curr_pred[9]).toFixed(2)}, Overlap: ${curr_pred[10]}<br>`
    curr_pred_str += `<a href="${gse_gene_viewer}"><button class="btn btn-primary btn-group-sm mt-2 mb-2">${gse} Gene Viewer</button></a>`
    
    curr_pred_str += "</div>"
    curr_pred_str += "<div class='p-3' style='background-color: rgb(247, 243, 248); border-radius: 1rem;'>"
    curr_pred_str += curr_pred[11].replaceAll('\n', '<br>')
    curr_pred_str += "</div>"

    selecter_div.innerHTML = selecter;
    hypothesis_div.innerHTML = curr_pred_str;
}

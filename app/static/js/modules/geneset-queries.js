// QUERY DIABETES RELATED SIGNATURES FOR GENESET
// [GeneSet]->[Signatures]
export async function geneset_signatures(geneset, resultid) {
    var inputvalue = geneset.split(',').join('\n');
    console.log(geneset)
    console.log(inputvalue)
    var desc = ''
    await $.ajax({
        url: "getdiabetesenrich",
        type: "POST",
        data: { genelist: inputvalue, description: desc }
    }).done(function (response) {

        const data = response['data']['Diabetes_Perturbations_GEO_2022'];

        if (!data) {
            $(`#${resultid}`).html("<p class='text-center m-5'> No data found </p>");
            return;
        }
        const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>"


        var tabletext = "<table id='table-enrichr' class='styled-table' style:'width=100%; vertical-align:top;'><thead><tr><th>Rank</th><th>Term name</th><th>P-value</th><th>Z-score</th><th>Combined score</th><th>Overlapping genes</th><th>Adjusted p-value</th></tr><tbody>";
        const currURL = window.location.href.split('/')
        for (var k = 0; k < data.length; k++) {
            tabletext += "<tr><td>" + data[k][0] + "</td><td>" + data[k][1] + "</td><td>" + Number(data[k][2]).toPrecision(4) + "</td><td>" + Number(data[k][3]).toPrecision(4) + "</td><td>" + Number(data[k][4]).toPrecision(4) + "</td><td>"
            var url = currURL.join('/') + 'singlegene'
            var api2 = currURL.join('/') + 'geneset'
            var gene_arr = data[k][5].map(g => `<a href='${url}' onclick="setGene('${g}')" target='_blank'>${g}<a/>`);

            tabletext += `<button class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-start text-left"
                    data-toggle="collapse" data-target="#genesoverlap-${data[k][0]}" aria-expanded="false"
                    aria-controls="genesoverlap-${data[k][0]}">
                    <div class="text">Show Overlapping Genes</div>
            </button>
                <div class="collapse" id="genesoverlap-${data[k][0]}">
                    ${gene_arr.join(", ")}
                    <a href="${api2}" onclick="setGenes('${data[k][5].join('&')}')" target='_blank'>
                        <button class="btn btn-primary btn-group-sm d-flex align-items-start text-center" style="font-size: small;">
                            Submit to Gene Set Queries
                        </button>
                    </a>
                </div></td>
            `
            tabletext += "<td>" + Number(data[k][6]).toPrecision(4) + "</td></tr>";
        }
        tabletext += "</tbody></table>";


        $(document).ready(function () {
            $('#table-enrichr').DataTable({
                dom: 'Bfrtip',
                buttons: [
                    'copy', { extend: 'csv', title: `${desc}-Diabetes-Perturbations-Enrichr-res` }
                ]
            });

        });
        document.getElementById(resultid).innerHTML = tabletext + clear_button;
    });
}

// QUERY DIABETES RELATED SIGNATURES FOR GENESET
// [GeneSet]->[Enrichment]
export async function geneset_enrichment(geneset, resultid) {
    var genes = geneset.split(',').join('\n')
    options = {}
    options.list = genes
    options.description = ''
    options.popup = true
    enrich(options)
}


// QUERY DIABETES RELATED SIGNATURES FOR GENESET
// [GeneSet]->[Kinases]
export async function geneset_kea3(geneset, resultid) {
    var inputvalue = geneset;
    await $.ajax({
        url: "getkea3",
        type: "POST",
        data: { "geneset": inputvalue}
    }).done(function (response) {
        console.log(response)
        const data = response["Integrated--meanRank"];

        if (!data) {
            $(`#${resultid}`).html("<p class='text-center m-5'> No data found </p>");
            return;
        }
        const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>"


        var tabletext = "<table id='table-kea3' class='styled-table' style:'width=100%; vertical-align:top;'><thead><tr><th>Rank</th><th>Kinase</th><th>Mean Rank</th><th>Overlapping genes</th></tr><tbody>";
        const currURL = window.location.href.split('/')
        for (var k = 0; k < data.length; k++) {
            tabletext += "<tr><td>" + data[k]['Rank'] + "</td><td>" + data[k]['TF'] + "</td><td>" + data[k]['Score'] + "</td><td>"
            var url = currURL.join('/') + 'singlegene'
            var api2 = currURL.join('/') + 'geneset'
            console.log(data[k].Overlapping_Genes)
            var gene_arr = data[k]['Overlapping_Genes'].split(',').map(g => `<a href='${url}' onclick="setGene('${g}')" target='_blank'>${g}<a/>`);

            tabletext += `<button class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-start text-left"
                    data-toggle="collapse" data-target="#genesoverlap-${data[k]['Rank']}" aria-expanded="false"
                    aria-controls="genesoverlap-${data[k]['Rank']}">
                    <div class="text">Show Overlapping Genes</div>
            </button>
                <div class="collapse" id="genesoverlap-${data[k]['Rank']}">
                    ${gene_arr.join(", ")}
                    <a href="${api2}" onclick="setGenes('${data[k]['Overlapping_Genes'].split(',').join('&')}')" target='_blank'>
                        <button class="btn btn-primary btn-group-sm d-flex align-items-start text-center" style="font-size: small;">
                            Submit to Gene Set Queries
                        </button>
                    </a>
                </div></td>
            `
        }
        tabletext += "</tbody></table>";


        $(document).ready(function () {
            $('#table-kea3').DataTable({
                dom: 'Bfrtip',
                buttons: [
                    'copy', { extend: 'csv', title: `kea3-results` }
                ]
            });

        });
        var appyter_url = `https://appyters.maayanlab.cloud/KEA3_Appyter/#/?args.Input%20gene/protein%20list=${inputvalue.split(',').join("%0A")}&submit`;
        var appyter_button = `
        <a id="kea3" target="_blank" rel="noopener noreferrer"
         href="${appyter_url}"><button type="button"
        class="btn btn-primary btn-group-sm mt-3 mb-3">
        <span id="appyter-action" class="ml-3">Start a new appyter in</span>
        <img src="static/img/appyters_logo.svg" class="img-fluid mr-3"
        style="width: 120px" alt="Appyters">
        </button>
        </a>`
        document.getElementById(resultid).innerHTML = tabletext + appyter_button + clear_button;
    });
}

// QUERY DIABETES RELATED SIGNATURES FOR GENESET
// [GeneSet]->[TFs]
export async function geneset_chea3(geneset, resultid) {
    var inputvalue = geneset;
    await $.ajax({
        url: "getchea3",
        type: "POST",
        data: { "geneset": inputvalue}
    }).done(function (response) {
        console.log(response)
        const data = response["Integrated--meanRank"];

        if (!data) {
            $(`#${resultid}`).html("<p class='text-center m-5'> No data found </p>");
            return;
        }
        const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>"


        var tabletext = "<table id='table-kea3' class='styled-table' style:'width=100%; vertical-align:top;'><thead><tr><th>Rank</th><th>Kinase</th><th>Mean Rank</th><th>Overlapping genes</th></tr><tbody>";
        const currURL = window.location.href.split('/')
        for (var k = 0; k < data.length; k++) {
            tabletext += "<tr><td>" + data[k]['Rank'] + "</td><td>" + data[k]['TF'] + "</td><td>" + data[k]['Score'] + "</td><td>"
            var url = currURL.join('/') + 'singlegene'
            var api2 = currURL.join('/') + 'geneset'
            console.log(data[k].Overlapping_Genes)
            var gene_arr = data[k]['Overlapping_Genes'].split(',').map(g => `<a href='${url}' onclick="setGene('${g}')" target='_blank'>${g}<a/>`);

            tabletext += `<button class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-start text-left"
                    data-toggle="collapse" data-target="#genesoverlap-${data[k]['Rank']}" aria-expanded="false"
                    aria-controls="genesoverlap-${data[k]['Rank']}">
                    <div class="text">Show Overlapping Genes</div>
            </button>
                <div class="collapse" id="genesoverlap-${data[k]['Rank']}">
                    ${gene_arr.join(", ")}
                    <a href="${api2}" onclick="setGenes('${data[k]['Overlapping_Genes'].split(',').join('&')}')" target='_blank'>
                        <button class="btn btn-primary btn-group-sm d-flex align-items-start text-center" style="font-size: small;">
                            Submit to Gene Set Queries
                        </button>
                    </a>
                </div></td>
            `
        }
        tabletext += "</tbody></table>";


        $(document).ready(function () {
            $('#table-kea3').DataTable({
                dom: 'Bfrtip',
                buttons: [
                    'copy', { extend: 'csv', title: `kea3-results` }
                ]
            });

        });
        var appyter_url = `https://appyters.maayanlab.cloud/ChEA3_Appyter/#/?args.paste_gene_input=${inputvalue.split(',').join("%0A")}&submit`;
        var appyter_button = `
        <a id="kea3" target="_blank" rel="noopener noreferrer"
         href="${appyter_url}"><button type="button"
        class="btn btn-primary btn-group-sm mt-3 mb-3">
        <span id="appyter-action" class="ml-3">Start a new appyter in</span>
        <img src="static/img/appyters_logo.svg" class="img-fluid mr-3"
        style="width: 120px" alt="Appyters">
        </button>
        </a>`
        document.getElementById(resultid).innerHTML = tabletext + appyter_button + clear_button;
    });
}

// QUERY DIABETES RELATED SIGNATURES FOR GENESET
// [GeneSet]->[SigComLincs]
export async function geneset_sigcomlincs(geneset, geneset_up, geneset_down, resultid) {
    var genes;
    if (geneset.length > 0) {genes = JSON.stringify({ 'genes': [geneset.split(',')] });}
    else {genes = JSON.stringify({ 'genes': [geneset_up.split(','), geneset_down.split(',')] });}
    console.log(genes)
    $.ajax({
        url: "getsigcom",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: genes
    }).done(function (response) {

        const url = response['url'];
        window.open(url, '_blank');

    });
}
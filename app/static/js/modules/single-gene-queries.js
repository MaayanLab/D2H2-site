export const human_list = fetch(
    "static/data/allgenes-comb.json"
).then(data => data.json());



export function gen_table(link, data, table_id, title, gene) {
    var titletext = `<div class ="row text-center mt-3"> <h4>${title}</h4></div>`;
    var tabletext = `<table id='${table_id}' class='styled-table' style='width: 90%;'><thead><tr>`
    tabletext += "<th>Signature</th><th>GEO Entry</th><th>P-value</th><th>Log2 Fold Change</th><th>Gene Rank in Signature</th><th>Boxplot Viewer</th></tr><tbody>"
    data.forEach(function (row) {
        var gse = row[row.length - 1].split("=")[1]
        var curr = window.location.href.replace('singlegene', '')
        var studyviewer = curr + gse
        tabletext += "<tr><td>" + row[0].split(/[-_]+/).join(" ") + "</td><td><a href='" + row[row.length - 1] + "' target='_blank'>" + gse + "</a></td><td>" + row[1] + "</td><td>" + row[2] + "</td><td>" + row[row.length - 2] + "</td><td><a href='" + studyviewer + `' target='_blank'><button class='btn btn-primary btn-group-sm' onclick="setGene('${gene}')">` + gse + " Gene Viewer</button></a></td></tr>"
    });
    tabletext += "</tbody></table>";
    var filename = link.split("/")[link.split("/").length - 1]
    document.getElementById("t2d-tables").innerHTML += (titletext + tabletext)
    $(document).ready(function () {
        $(`#${table_id}`).DataTable({
            order: [[2, 'asc']],
            dom: 'Bfrtip',
            buttons: [
                'copy', { extend: 'csv', title: filename }
            ]
        });
    })
}


// SINGLE GENE DIABETES SIGNATURES:
// [Gene]->[Signatures]
export async function gene_signatures(gene, species, resultid) {
    var clear_button;
    if (resultid.includes('[')) clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-4 mb-1' onclick='clear_home();'> Clear Results </button> </a>";
    else clear_button = "";

    var jsonData = {};
    species = species.toLowerCase();
    jsonData["gene"] = gene;
    jsonData["species"] = species;
    document.getElementById(resultid).innerHTML = `<div class='container justify-content-center'><div id='volcano-plot-${resultid}' class='justify-content-center centered'></div><div id='buttons'></div><div id='t2d-tables'></div></div>`
    await $.ajax({
        url: "api/volcano",
        type: "POST",
        data: jsonData,
    }).then(function (res) {
        var jdata = JSON.parse(res)
        var plot = jdata['plot']
        var tables = jdata['tables']
        var table_values = jdata['table_values']
        var micro = jdata['micro']
        plot['target_id']= `volcano-plot-${resultid}`;
        console.log(plot)
        if (document.getElementById("buttons")) {
            document.getElementById("buttons").innerHTML += `<div class='row text-center justify-content-center'>${clear_button}</div>`
        }
        window.Bokeh.embed.embed_item(plot)
        var dir = "up";
        var titleRNA = `Top ${species} RNA-seq signatures where ${gene} is ${dir}-regulated`
        var titlemicro = `Top ${species} microarray signatures where ${gene} is ${dir}-regulated`

        gen_table(tables[0], table_values[0], `${species}_up`, titleRNA, gene)
        if (micro) gen_table(tables[2], table_values[2], `${species}_micro_up`, titlemicro, gene)
        dir = "down";
        var titleRNA = `Top ${species} RNA-seq signatures where ${gene} is ${dir}-regulated`
        var titlemicro = `Top ${species} microarray signatures where ${gene} is ${dir}-regulated`
        gen_table(tables[1], table_values[1], `${species}_down`, titleRNA, gene)
        if (micro) gen_table(tables[3], table_values[3], `${species}_micro_down`, titlemicro, gene)
    });
}


// SINGLE GENE EXPRESSION:
// [Gene]->[Expression]
export async function generanger_plot(gene, library, resultid) {
    document.getElementById(resultid).innerHTML = `<div class='container justify-content-center mb-5 mx-auto text-center' style='overflow: scroll;'><div id='volcano-plot-${resultid}' class='justify-content-center centered'></div><div id='buttons' class='text-center'></div></div>`
    var jsonData = {};
    jsonData["gene"] = gene;
    await $.ajax({
        url: "queryexpression",
        type: "POST",
        dataType: 'json',
        data: jsonData,
        success: function (jdata) {
            if (Object.keys(jdata.allData.dbData).length == 0) {
                document.getElementById(resultid).innerHTML = "<p class='text-center'>No data for this gene found</p>"
                return;
            }
            var plotLib = jdata['allData']['dbData'][library]
            if (plotLib) {
            const val_names = ['lowerfence', 'upperfence', 'mean', 'median', 'q1', 'q3', 'sd', 'names', 'levels']
            val_names.forEach((attr) => {
                if (attr in plotLib) plotLib[attr] = plotLib[attr].reverse();
            });
            if ('levels' in plotLib) {
                plotLib['type'] = 'scatter'
                plotLib['marker'] = {color: '#1f77b4'}
                plotLib["mode"]= 'markers'
                plotLib["y"]= plotLib["levels"]
                plotLib["x"]= plotLib["names"]
            } else {
                plotLib.x = plotLib.names
                plotLib.type = 'box'
                plotLib.orientation = 'v'
            }
            let customWidth = plotLib.x.length * 20;
            var layout = {
                autosize: true,
                width: customWidth,
                height: 600,
                title: {text:  gene + ' '+ library + ' Expression', xanchor: 'left', yanchor: 'top', x: 0},
                xaxis: {
                    automargin: true,
                    tickangle: 45,
                    range: [-0.5, plotLib['x'].length]
                }
            };
            Plotly.newPlot(`volcano-plot-${resultid}`, [plotLib], layout)

            document.getElementById(resultid).innerHTML += `
            <div='row'>
            <a id="gtex-url" target="_blank" rel="noopener noreferrer" href=https://generanger.maayanlab.cloud/gene/${gene}?database=${library}>
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                <img src="static/img/generangerlogo.png" class="img-fluid ml-1 mr-3" style="width: 30px" alt="GTEx">
                </button>
            </a>
            <a id="archs-url" target="_blank" rel="noopener noreferrer" href="https://maayanlab.cloud/archs4/gene/${gene}">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                <img src="static/img/archs4logo.png" class="img-fluid mr-3" style="width: 110px" alt="ARCHS4">
                </button>
            </a>
            <a id="gtex-url" target="_blank" rel="noopener noreferrer" href="https://gtexportal.org/home/gene/${gene}">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                <img src="static/img/gtexlogo.png" class="img-fluid mr-3" style="width: 110px" alt="GTEx">
                </button>
            </a>
            </div>
            `
        } else {
            document.getElementById(resultid).innerHTML = "<p>I'm sorry, no data for this gene was available from the chosen library.</p>"
        }
        }
    });
}


// SINGLE GENE PERTURBATIONS:
// [Gene]->[Perturbations]
export async function geo_reverse(gene, species) {
        if (species == 'Human') {
            var arg = 'human_gene';
        } else {
            var arg = 'mouse_gene';
        }
        const formData = new FormData()
        formData.append('species_input', species)
        formData.append(arg, gene)
        var res = await fetch("https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/", {
            method: "POST",
            headers: {
                'Accept': 'application/json',
            },
            body: formData,
        })
        const id = await res.json()
        window.open("https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/" + id.session_id, '_blank')
}

export async function l1000_reverse(gene) {
    var check_list = await human_list
    if (!check_list.includes(gene)) {
        alert('This Appyter only accepts Human gene symbols.')
        return;
    }
    window.open("https://appyters.maayanlab.cloud/L1000_RNAseq_Gene_Search/#/?args.gene=" + gene + "&submit", '_blank');
}

export function single_gene_perturbations(gene, species, resultid) {
    document.getElementById(resultid).innerHTML = `
              <div class="row" style="flex-wrap: wrap;">
                <div class="col-md-6 col-lg-6 col-sm-12">
                  <div class="row text-center">
                    <p class="mt-4 ml-1">This Appyter can be used to find conditions to maximally up/down regulate the expression of a gene in human/mouse based on curated GEO studies.</p>
                  </div>
                  <div class="row justify-content-center mt-3 mb-3">
                    <img src="static/img/genecentric-thumbnail.png" alt="" class="img-fluid img-thumbnail">
                  </div>
                  <div class="justify-content-center row text-center mb-3">
                    <a id="geo_reverse" target="_blank" rel="noopener noreferrer"><button type="button"
                        class="btn btn-primary btn-group-sm mt-3 mb-3">
                        <span id="appyter-action" class="ml-3">Start a new appyter in</span>
                        <img src="static/img/appyters_logo.svg" class="img-fluid mr-3" style="width: 120px"
                          alt="Appyters">
                      </button>
                    </a>
                  </div>
                </div>
                <div class="col-md-6 col-lg-6 col-sm-12">
                  <div class="row text-center mr-2">
                    <p class="mt-3">This Appyter generates a volcano plot displaying how different drugs and small molecules may induce or suppress the expression of a specific gene based on transformed L1000 data.</p>
                  </div>
                  <div class="row justify-content-center mb-3">
                    <img src="static/img/L1000-thumbnail.png" alt="" class="img-fluid img-thumbnail">
                  </div>
                  <div class="justify-content-center row text-center mb-3">
                    <a id="l1000_reverse" target="_blank" rel="noopener noreferrer"><button type="button"
                        class="btn btn-primary btn-group-sm mt-3 mb-3">
                        <span id="appyter-action" class="ml-3">Start a new appyter in</span>
                        <img src="static/img/appyters_logo.svg" class="img-fluid mr-3" style="width: 120px"
                          alt="Appyters">
                      </button>
                    </a>
                  </div>
                </div>
              </div>`
    document.getElementById("l1000_reverse").addEventListener("click",  function() {
        l1000_reverse(gene);
    });
    document.getElementById("geo_reverse").addEventListener("click",  function() {
        geo_reverse(gene, species);
    });
}

// SINGLE GENE TRAITS:
// [Gene]->[Traits]
export async function query_gwas(gene, id) {
    var inputvalue = gene;
    await $.ajax({
        url: "getgwas",
        type: "POST",
        data: { gene: inputvalue }
    }).done(function (response) {
        const data = response['GWAS_Catalog'];
        if (data.length === 0) {
            $(`#${id}`).html("<p class='text-center'> No data found </p>");
            return;
        }
        var clear_button;
        if (id.includes('[')) {
            clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>";
        } else {
            clear_button = "";
        }
        
        var tabletext = `<table id='table-gwas${id}' class='styled-table'><thead><tr><th></th><th>Gene</th><th>Trait</th><th>Count</th></tr><tbody>`;
        for (var k = 0; k < data.length; k++) {
            tabletext += "<tr><td>" + (k + 1) + "</td><td><a href='https://www.ebi.ac.uk/gwas/genes/" + data[k]['gene'] + "' target='_blank'>" + data[k]['gene'] + "</a></td><td><a href='" + data[k]['mapped_trait_link'] + "' target='_blank'>" + data[k]['trait'] + "</a></td><td>" + data[k]['count'] + "</td></tr>";
        }
        tabletext += "</tbody></table>";
        $(document).ready(function () {
            $(`#table-gwas${id}`).DataTable({
                dom: 'Bfrtip',
                buttons: [
                    'copy', { extend: 'csv', title: `${inputvalue}-gwas-res` }
                ]
            });
        });
        document.getElementById(id).innerHTML = tabletext + clear_button;
    });
}

// SINGLE GENE TFS
// [Gene]->[TFs]
export async function query_enrichr_tfs(gene, id) {
    var inputvalue = gene;
    if (!inputvalue) {
        document.getElementById(id).innerHTML = "<p class='text-center'> No data found </p>";
        return;
    }
    await $.ajax({
        url: "gettfs",
        type: "POST",
        data: { gene: inputvalue }
    }).done(function (response) {
        const data = response['data'];
        var clear_button;
        if (id.includes('[')) {
            clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>";
        } else clear_button = "";
        
        if (data.length === 0) {
            document.getElementById(id).innerHTML = "<p class='text-center'> No data found </p>";
            return;
        }
        var selecter = `<div class='text-center'><p>Select from one the annotated libraries: </p><select class="m-2 libpicker" data-style="btn-primary" onchange="on_change(this)" data-width="500px">`
        for (var i = 0; i < data.length; i++) {
            var lib = data[i]['name'];
            var libDisplay = lib.replaceAll("_", " ")
            selecter += `<option value=${lib}>${libDisplay}</option>`;
        }
        selecter += `</select></div>`
        for (var i = 0; i < data.length; i++) {
            var lib = data[i]['name'];
            var sentence = data[i]['format'];
            var res_html = `<div id="${lib}" style="display:none;"><table id='table-enrichr' class='styled-table table-enrichr'><thead><tr><th></th></tr><tbody>`
            for (var j = 0; j < data[i]['tfs'].length; j++) {
                var tf = data[i]['tfs'][j]
                var tf_sentence = sentence.replace('{0}', inputvalue).replace('{1}', tf)
                res_html += `<tr><td> ${tf_sentence} </td></tr>`
            }
            res_html += `</tbody></table></div>`
            selecter += res_html
        }
        selecter += `<script>      
                    </script>`
        $(document).ready(function () {
            $(`.table-enrichr`).DataTable({
                dom: 'Bfrtip',
                buttons: [
                    'copy', { extend: 'csv', title: `${inputvalue}-enrichr-tfs` }
                ]
            });
        });
        document.getElementById(id).innerHTML = selecter + clear_button;
        document.getElementById(data[0]['name']).style.display = 'block';
    });
}

// SINGLE GENE QUERY ARCHS4 FOR CORRELATIONS
// [Gene]->[Correlation]
export async function loadCorrelation(gene, id) {
    var jsonData = {};

    jsonData["id"] = gene;
    jsonData["count"] = 101;
    await $.ajax({
        type: "POST",
        url: "https://maayanlab.cloud/matrixapi/coltop",
        contentType: "application/json; charset=utf-8",
        dataType: "json",
        data: JSON.stringify(jsonData),
        success: function (jdata) {
            var data = jdata;

            var genesym = data["rowids"];
            var correlation = data["values"];

            if (!(genesym)) {
                $(`#${id}`).html("<p class='text-center'> No data found </p>");
                return;
            }
            var clear_button;
            if (id.includes('[')) {
                clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>";
            } else clear_button = "";

            //3) separate them back out:
            var genes = "";
            var tabletext = "<table id='tablecor' class='styled-table'><thead><tr><th>Rank</th><th>Gene Symbol</th><th>Pearson Correlation</th></tr><tbody>";
            for (var k = 1; k < genesym.length; k++) {
                tabletext += "<tr><td>" + k + "</td><td><a href=\"https://maayanlab.cloud/archs4/gene/" + genesym[k] + "\" target=\"_blank\">" + genesym[k] + "</a></td><td>" + Number(correlation[k]).toPrecision(4) + "</td></tr>";
                genes = genes + genesym[k] + "\n";
            }
            tabletext += "</tbody></table>";

            $(document).ready(function () {
                $('#tablecor').DataTable({
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: `${gene}-archs4-corr` }
                    ]
                });
            });

            document.getElementById(id).innerHTML = tabletext + clear_button;

        },
        error: function (xhr, textStatus, errorThrown) {
        }
    });
}

// SINGLE GENE QUERY ARCHS4 FOR CORRELATIONS
// [Gene]->[Knockout]
export async function query_komp(gene, id) {
    var inputvalue = gene;
    await $.ajax({
        url: "getkomp",
        type: "POST",
        data: { gene: inputvalue }
    }).done(function (response) {
        const data = response['data'];
        if (data.length === 0) {
            $(`#${id}`).html("<p class='text-center'> No data found </p>");
            return;
        }
        
        var clear_button;
        if (id.includes('[')) {
            clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_home();'> Clear Results </button> </a>";
        } else clear_button = '';
        var tabletext = "<table id='table-pheno' class='styled-table'><thead><tr><th>Gene</th><th>Phenotype</th><th>PM ID</th><th>Comments</th></tr><tbody>";
        for (var k = 0; k < data.length; k++) {
            tabletext += "<tr><td><a href='http://www.informatics.jax.org/marker/" + data[k]['OntologyAnnotation.subject.primaryIdentifier'] + "' target='_blank'>" + data[k]['OntologyAnnotation.subject.symbol'] + "</a></td>"
            tabletext += "<td><a href='http://www.informatics.jax.org/vocab/mp_ontology/" + data[k]['OntologyAnnotation.ontologyTerm.identifier'] + "' target='_blank'>" + data[k]['OntologyAnnotation.ontologyTerm.name'] + "</a></td>"
            tabletext += "<td><a href='https://pubmed.ncbi.nlm.nih.gov/" + data[k]['OntologyAnnotation.evidence.publications.pubMedId'] + "' target='_blank'>" + data[k]['OntologyAnnotation.evidence.publications.pubMedId'] + "</a></td>"
            tabletext += "<td>" + data[k]['OntologyAnnotation.evidence.comments.description'] + "</td></tr>"
        }
        tabletext += "</tbody></table>";

        $(document).ready(function () {
            $('#table-pheno').DataTable({
                dom: 'Bfrtip',
                buttons: [
                    'copy', { extend: 'csv', title: `${inputvalue}-mgi-pheno` }
                ]
            });
        });
        document.getElementById(id).innerHTML = tabletext + clear_button;
    });
}
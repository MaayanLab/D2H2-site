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

export async function search_for_genesets(term, resource, id) {
    if (resource == 'Enrichr') {
        genesets = query_enrichr_metadata(term)
    } else if (resource == 'Rummagene') {
        query_rummagene(term, id)

    } else if (resource == 'Geneshot') {

    }

    return 
}


// Term to Enrichr Metadata
// [Term]->[Enrichr]
export async function query_enrichr_metadata(term) {
    const {terms} = await (
        await fetch(`https://maayanlab.cloud/Enrichr/termmap?meta=${term}&json=true`, {
            method: 'GET'
        })
    ).json()
    const options = []
    for (const [library,v] of Object.entries(terms)) {
        for (const term of v)  {
            options.push({
                library,
                term
            })
        }
    }
    console.log(options)

    

    return options
}


export async function get_enrichr_geneset(term, library) {
    const res = await fetch(`https://maayanlab.cloud/Enrichr/geneSetLibrary?term=${term}&libraryName=${library}&mode=json`, {
        method: 'GET',
    })
    const results = await res.json()
    const vals = Object.entries(results)           
    if (vals) {
        const [description, genes] = vals[0]
        return {description: `${description} (${library})`, genes}
    }
}

export async function query_rummagene(term) {

    /* const query = `query MyQuery {
        geneSetTermSearch(terms: "${term}", first: 1000) {
          edges {
            node {
              term
              id
              nGeneIds
            }
          }
        }
      }`

    const res = await fetch("https://rummagene.com/graphiql", {
        method: "POST",
        body: query,
        headers: {
            "Content-type": "application/graphiql"
        }
    }).json()

    console.log(res) */

    const url1 = `https://rummagene.com/pubmed-search?q=${term}`;
    const url2 = `https://rummagene.com/term-search?q=${term})`;


    document.getElementById(resultid).innerHTML =
        `
        <div class="row">
            Search by querying column names:
            <a href="${url1}" target="_blank">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                    <img src="static/img/rummagene_logo.png" class="img-fluid mr-3"
                        style="width: 60px" alt="Rummagene">
                </button>
            </a>
            Search by querying PMC:
            <a href="${url2}" target="_blank">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                    <img src="static/img/rummagene_logo.png" class="img-fluid mr-3"
                        style="width: 60px" alt="Rummagene">
                </button>
            </a>
        </div>
        `
}

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


// Term to Enrichr Metadata
// [Term]->[Enrichr]
export async function query_enrichr_metadata(term, id) {
    const {terms} = await (
        await fetch(`https://maayanlab.cloud/Enrichr/termmap?meta=${term}&json=true`, {
            method: 'GET',
            signal: controller.signal
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
}







const {library, term} = params
const controller = get_controller() 
const res = await fetch(`https://maayanlab.cloud/Enrichr/geneSetLibrary?term=${term}&libraryName=${library}&mode=json`, {
    method: 'GET',
})
const results = await res.json()
const vals = Object.entries(results)           
if (vals) {
    const [description, genes] = vals[0]
    setInput({description: `${description} (${library})`, genes})
    setOpen(false)
}


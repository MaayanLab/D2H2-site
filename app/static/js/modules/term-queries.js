export async function search_for_genesets(term, resource, id) {
    if (resource == 'Enrichr') {
        await query_enrichr_metadata(term, id);
    } else if (resource == 'Rummagene') {
        await query_rummagene(term, id);
    } else if (resource == 'Geneshot') {
        await query_geneshot(term, id);
    }
    return;
}


// Term to Enrichr Metadata
// [Term]->[Enrichr]
export async function query_enrichr_metadata(term, id) {
    const {terms} = await (
        await fetch(`https://maayanlab.cloud/Enrichr/termmap?meta=${term}&json=true`, {
            method: 'GET'
        })
    ).json()
    const options = []
    var libraries = new Set();
    for (const [library,v] of Object.entries(terms)) {
        for (const t of v)  {
            if (t.toLowerCase().includes(term.toLowerCase())) {
                libraries.add(library)
                options.push({
                    library,
                    t
                })
            }
        }
    }

    if (options.length < 1) {
        document.getElementById(id).innerHTML = "No gene sets were found from Enrichr containing your term.";
        return
    }

    libraries = Array.from(libraries)

    var selecter = `<div class='text-center'><p>Select from one the annotated libraries: </p><select class="m-2 libpicker" data-style="btn-primary" onchange="on_change(this)" data-width="500px">`
    for (var i = 0; i < libraries.length; i++) {
        var lib = libraries[i];
        var libDisplay = lib.replaceAll("_", " ")
        selecter += `<option value=${lib}>${libDisplay}</option>`;
    }
    selecter += `</select></div>`
    for (var i = 0; i < options.length; i++) {
        var lib = options[i].library;
        const lib_terms = options.filter(l => l.library == lib)
        var res_html = `<div id="${lib}" style="display:none;"><table id='table-enrichr' class='styled-table table-enrichr'><thead><tr><th></th><th></th></tr><tbody>`
        for (var j = 0; j < lib_terms.length; j++) {
            var t = lib_terms[j].t
            res_html += `<tr><td> ${t} </td><td> <button type="button" class="btn btn-primary btn-group-sm" onclick="show_geneset_modal('${lib}', '${t}')">View Geneset</button> </td></tr>`
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
                'copy', { extend: 'csv', title: `${term}-enrichr-genesets` }
            ]
        });
    });

    const placeholder = document.createElement("div");
    placeholder.innerHTML = selecter;
    document.getElementById(id).appendChild(placeholder);
    document.getElementById(options[0].library).style.display = 'block';
    return; 
}


export async function query_rummagene(term, id) {

    const url1 = `https://rummagene.com/pubmed-search?q=${term}`;
    const url2 = `https://rummagene.com/term-search?q=${term}`;


    document.getElementById(id).innerHTML =
        `
        <div class="col text-right">
            <div class="col">
            Search by querying PMC:
            <a href="${url1}" target="_blank">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                    <img src="static/img/rummagene_logo.png" class="img-fluid mr-3"
                        style="width: 60px" alt="Rummagene">
                </button>
            </a>
            </div>
            <div class="col">
            Search by querying column names:
            <a href="${url2}" target="_blank">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"> Open in
                    <img src="static/img/rummagene_logo.png" class="img-fluid mr-3"
                        style="width: 60px" alt="Rummagene">
                </button>
            </a>
            </div>
        </div>
        `
}

export async function query_geneshot(term, id) {
    var data = JSON.stringify({'term': term })
    const res = await $.ajax({
        url: "api/run_geneshot",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: data
    })

    if (res.hasOwnProperty("error")) {
        document.getElementById(id).innerHTML = '<div>The Geneshot API is currently down. Please try again later.</div>';
        return;
    }

    const file = new Blob([`${res.description}\t${res['genes'].join('\t')}`], { type: 'text/plain' });
    // Add file content in the object URL
    const href = URL.createObjectURL(file);
    const download = `${term}_Geneshot.txt`;

    const result = `
    <div id="geneshot-text-area" style="margin: auto;">
        <div class="text-center justify-content-center">
            <p id="geneshot-title">Geneshot produced a gene set based on literature comentions with following term ${term} (${res['count']})</p>
            <textarea id="geneshot-text" class="text-center" style="height: 20rem; overflow-y: scroll; white-space: pre-wrap;" readonly>
            ${res['genes'].join('\r\n')}
            </textarea>
        </div>
        <div class="row">
            <button type="button" onclick="navigator.clipboard.writeText(document.getElementById('geneshot-text').value)" class="btn btn-primary btn-group-sm m-2 ">
                Copy to clipboard
            </button>
            <a id="geneshot-queries" href="{{ url_for('geneset_home') }}" onclick="setGenes('${res['genes'].join('&')}')" target='_blank' style="text-decoration: none">
                <button type="button" class="btn btn-primary btn-group-sm m-2">
                    Open in Geneset Queries
                </button>
            </a>
            <a id="geneshot-download" href="${href}" download="${download}" style="text-decoration: none">
                <button type="button" class="btn btn-primary btn-group-sm m-2">
                Download
                </button>
            </a>
        </div>
    </div>`
    document.getElementById(id).innerHTML = result;
}


function check_genes_present(genes) {
    if (genes.length == 0) {
        alert('no significantly differentially expressed genes were identified with these thresholds');
        return true;
    } else return false;
}

function filter_genes_sigs(genelist, adjpvals, pvals) {
    var col = 'pval'
    try {
        col = document.getElementById('col-to-use').value
        console.log(col)
    } catch {
        col = 'pval'
    }

    genelist = genelist.split(',')
    if (col == 'pval'){
        var sigs = pvals;

    } else var sigs = adjpvals;
    sigs = sigs.split(',').map(function(item) {
        return parseFloat(item);
    });
    var numgenes = document.getElementById('numgenes').value
    var signifigance = document.getElementById('signifigance').value
    var genes_valid = []

    for (i=0; i < genelist.length; i++ ){
        if (sigs[i] <= signifigance) {
            genes_valid.push(genelist[i])
        }
    }
    console.log(genes_valid)
    return genes_valid
}


function generate_single_plots(){
    // This function will generate the umap, tsne, and pca plots for each indivdual study for a specific condition
    document.getElementById("umap-plot").innerHTML = "";
    document.getElementById("tsne-plot").innerHTML = "";
    document.getElementById("pca-plot").innerHTML = "";

    

    // document.getElementById("singleplots-loading").innerHTML = "<div class='loader justify-content-center'></div>";
    $('#singleplots-loading').addClass('loader justify-content-center');
    var gse = document.getElementById("singlegse").innerText
    var species = document.getElementById("species").innerText
    var condition_group = document.getElementById("methodsingle").value

    var gsedata = JSON.stringify({'gse': gse, 'species': species, 'conditiongroup':condition_group});
    $.ajax({
        url: "/singleplots",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: gsedata,
    }).done(function(response) {
        document.getElementById("singleplots-loading").innerHTML = "";

        $('#singleplots-loading').removeClass('loader justify-content-center');
        var umap_plot = response['umapplot']
        var pca_plot = response['pcaplot']
        var tsne_plot = response['tsneplot']

        window.Bokeh.embed.embed_item(umap_plot)
        window.Bokeh.embed.embed_item(pca_plot)
        window.Bokeh.embed.embed_item(tsne_plot)
    });
    $("#boxplot").attr("data-url-plot", `/api/plot_single/${gse}/${condition_group}`)
}


const human_list = fetch(
    "static/searchdata/allgenes-comb.json"
).then(data => data.json());                 


///////// anitmated number counters /////////////
// How long you want the animation to take, in ms
const animationDuration = 2000;
// Calculate how long each ‘frame’ should last if we want to update the animation 60 times per second
const frameDuration = 1000 / 60;
// Use that to calculate how many frames we need to complete the animation
const totalFrames = Math.round( animationDuration / frameDuration );
// An ease-out function that slows the count as it progresses
const easeOutQuad = t => t * ( 2 - t );

// The animation function, which takes an Element
const animateCountUp = el => {
  let frame = 0;
  const countTo = parseInt( el.innerHTML, 10 );
  // Start the animation running 60 times per second
  const counter = setInterval( () => {
    frame++;
    // Calculate our progress as a value between 0 and 1
    // Pass that value to our easing function to get our
    // progress on a curve
    const progress = easeOutQuad( frame / totalFrames );
    // Use the progress value to calculate the current count
    const currentCount = Math.round( countTo * progress );

    // If the current count has changed, update the element
    if ( parseInt( el.innerHTML, 10 ) !== currentCount ) {
      el.innerHTML = currentCount;
    }

    // If we’ve reached our last frame, stop the animation
    if ( frame === totalFrames ) {
      clearInterval( counter );
    }
  }, frameDuration );
};


// COUNT GENES IN TEXT BOXS ON geneset page

function geneCount(gene_list, num) {
    const genes = gene_list.toUpperCase().split(/\r?\n/g).filter(Boolean);
    $('span#gene-count' + String(num)).text(genes.length);
}

function clear_home() {
    document.getElementById("t2d-tables").innerHTML = ""
    document.getElementById("volcano-plot").innerHTML = ""
    document.getElementById("buttons").innerHTML = ""

}

function submit_geneset(genelist, adjpvals, pvals) {

    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value
    
    
    if (check_genes_present(genes_valid)) return;


    var genes = genes_valid.splice(0, numgenes).join('&')
        



    const currURL = window.location.href.split('/')
    localStorage.setItem('geneset', genes)
    var url = currURL.filter(x => !x.includes('GSE')).join('/') + '/geneset'
    window.open(url, '_blank')
}

function submit_geneset_home(genelist, adjpvals, pvals, descset) {

    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value

    if (check_genes_present(genes_valid)) return;

    var genes = genes_valid.splice(0, numgenes).join('&')


    localStorage.setItem('genes', genes)
    localStorage.setItem('descset', '')
    var home = window.location.href.split('/').filter(x => !x.includes('GSE')).join('/')
    window.open(home, '_blank')
}



function setGenes(genes) {
    localStorage.setItem('geneset', genes)
}

function clear_dge() {
    document.getElementById("dge-table-area").innerHTML = ""
    document.getElementById("dge-plot").innerHTML = ""
    document.getElementById("dge-loading").innerHTML = ""
    document.getElementById("geneset-buttons").innerHTML = ""
    
}
function clear_dge_single() {
    document.getElementById("dge-table-area").innerHTML = ""
    document.getElementById("dge-plot").innerHTML = ""
    document.getElementById("dge-loading").innerHTML = ""
    document.getElementById("geneset-buttons").innerHTML = ""
    document.getElementById("enrichment-area").innerHTML = ""
    
}

function on_change(el) {

    for (var i =0; i < el.options.length; i++) {
        document.getElementById(el.options[i].value).style.display = 'none';
    }
    document.getElementById(el.options[el.selectedIndex].value).style.display = 'block'; // Show el

}


function gen_table(link, table_id, title, gene) {
    var csvdata = parseCsv(link)
    csvdata.then(function(data) {
        var titletext = `<div class ="row text-center mt-3"> <h4>${title}</h4></div>`;
        var tabletext = `<table id='${table_id}' class='styled-table'><thead><tr>`


        tabletext += "<th>Signature</th><th>GEO Entry</th><th>P-value</th><th>Log2 Fold Change</th><th>Gene Rank in Signature</th><th>Boxplot Viewer</th></tr><tbody>"


        data.data.forEach(function(row) {
            var gse = row['Link to GEO Study'].split("=")[1]
            var curr = window.location.href
            var studyviewer = curr + gse
            tabletext += "<tr><td>" +row['Signature'] +"</td><td><a href='" + row['Link to GEO Study'] + "' target='_blank'>" + gse +"</a></td><td>" +row['P-value'] +"</td><td>"+row['Log2 Fold Change'] + "</td><td>"+row['Gene Rank in Signature'] + "</td><td><a href='"+ studyviewer + `' target='_blank'><button class='btn btn-primary btn-group-sm' onclick="setGene('${gene}')">`+gse + " Gene Viewer</button></a></td></tr>"
        });

        tabletext += "</tbody></table>";
        var filename = link.split("/")[link.split("/").length -1]
        var download = `Download table: <a href="${link}">${filename}</a>`
        document.getElementById("t2d-tables").innerHTML += (titletext + tabletext + download)

        $(document).ready(function() {
            $(`#${table_id}`).DataTable({
                order: [[2, 'asc']],
            });
        })
        
    })
}

async function parseCsv(file) {
    return new Promise((resolve, reject) => {
        Papa.parse(file, {
            download: true,
            header: true,
            skipEmptyLines: true,
            complete: (results) => {
            return resolve(results);
            },
            error: (error ) => {
            return reject(error);
            },
        });
    });
}


// OPEN GENE LIST IN ENRICHR
function enrich(options) {
    if (typeof options.list === 'undefined') {
        alert('No genes defined.');
        return;
    }

    var description  = options.description || "",
    	popup = options.popup || false,
    	form = document.createElement('form'),
    	listField = document.createElement('input'),
    	descField = document.createElement('input');
  
    form.setAttribute('method', 'post');
    form.setAttribute('action', 'https://maayanlab.cloud/Enrichr/enrich');
    if (popup) {
        form.setAttribute('target', '_blank');
    }
    form.setAttribute('enctype', 'multipart/form-data');

    listField.setAttribute('type', 'hidden');
    listField.setAttribute('name', 'list');
    listField.setAttribute('value', options.list);
    form.appendChild(listField);

    descField.setAttribute('type', 'hidden');
    descField.setAttribute('name', 'description');
    descField.setAttribute('value', description);
    form.appendChild(descField);

    document.body.appendChild(form);
    form.submit();
    document.body.removeChild(form);
}

function filter_and_submit_to_enrichr(genelist, adjpvals, pvals, description) {
    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value
    var genes = genes_valid.splice(0, numgenes).join('\n')
        
    options = {}
    options.list = genes
    options.description = description
    options.popup = true
    enrich(options)
}

async function filter_and_submit_to_kg(genelist, adjpvals, pvals, description) {

    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value

    var genes_str = genes_valid.splice(0, numgenes).join('\n')

    const formData = new FormData()
    formData.append('list', genes_str)
    formData.append('description', description)

    var res = await fetch("https://maayanlab.cloud/Enrichr/addList", {
        method: "POST",
        headers: {
            'Accept': 'application/json',
        },
        body: formData,
    })

    const result = await res.json()
    console.log(result)
    const userListId = result['userListId']
    const url = `https://enrichr-kg.dev.maayanlab.cloud/?userListId=${userListId}`
    window.open(url, '_blank')
}


function loadFileAsText(section, delim) {
    
    return new Promise((resolve, reject) => {
    var fileToLoad = document.getElementById(section).files[0];
  
    var fileReader = new FileReader();
    fileReader.onload = async function(fileLoadedEvent){
        var textFromFileLoaded = fileLoadedEvent.target.result;
        var genes = textFromFileLoaded.split(/[\n\t,]/).join(delim)
        resolve(genes)
    };

    fileReader.readAsBinaryString(fileToLoad, "UTF-8");
    });   
}

function fillSingleExample(gene) {
    $(document).ready(function() {
        document.getElementById("singlegenenav").classList.add('active')
        for (var i = 1; i < 7; i++) {
            
            var selectize = $(`#search${i}`)[0].selectize;
            selectize.setValue(gene);
        }
    })
}

function fillSingleExampleSkip(gene, skip) {
    for (var i = 1; i < 7; i++) {
        if (`#search${i}` != skip) {
        var selectize = $(`#search${i}`)[0].selectize;
        selectize.setValue(gene);
        }
    }
}


function fillSingleExampleHome(gene) {
    var selectize = $(`#gene-select`)[0].selectize;
    selectize.setValue(gene);
}

function fillSetExample(geneset) {
    $('.input-form').each(function() {
        document.getElementById(this.id).value = geneset;
        geneCount(geneset, this.id[this.id.length -1])
    })
}


function fillSet(id, descid, count_id) {
    $.ajax({
        url: "getexample",
        type: "POST",
        data: {},
        dataType: 'json',
    }).done(function(response) {

        const desc = response['description']
        const genes = response['genes']
        document.getElementById(id).value = genes;
        if (descid != '') {
            document.getElementById(descid).value = desc;
        }
        geneCount(genes, count_id)

    });
}

function setGene(gene) {
    localStorage.setItem('gene', gene)
}


$(document).ready(function() {

    // SMALL NAV MENU

    $("<select class='selectize-nav text-center justify-content-center hamburger'> <img href='static/img/hamburger.jpeg/> </select>").appendTo("#mainnav");

    // Create default option "Go to..."
    $("<option />", {
       "selected": "selected",
       "value"   : "",
       "text"    : "Go to.."
    }).appendTo("nav select");
    
    // Populate dropdown with menu items
    $("nav a").each(function() {
     var el = $(this);
     if (el.attr("href").substr(0, 1) != "#" && el.attr("href") != '/' && el.attr("id") != 'toc') {
        $("<option />", {
            "value"   : el.attr("href"),
            "text"    : el.text().trim()
        }).appendTo("#mainnav select");
    }
    });

    $("nav select").change(function() {
        window.location = $(this).find("option:selected").val();
    });

    // BOLD CURRENT PAGE

    $('.nav-link').each(function(){
        var url = window.location.href
        if (url.includes('#')) {
            url = url.split('#')[0]
        }
    
        if ($(this).prop('href') == url) {
            $(this).addClass('active'); 
            $(this).parents('li').addClass('active');
        }
        // if (url.split('/')[3].startsWith('GSE')) {
        //     $("#viewer").addClass('active'); 
        //     $("#viewer").parents('li').addClass('active');
        // }
    });

    var currURL = window.location.href.split("/");


    var first_load = true;

    // SELECTIZE FOR VIEWER AND HOME PAGE APPYTER T2D

    var $gene_select = $('#gene-select').selectize({
        preload: true,
        valueField: 'gene_symbol',
        labelField: 'gene_symbol',
        searchField: 'gene_symbol',
        render: {
            option: function (item, escape) {
                return '<div class="pt-2 light">'+item.gene_symbol+'</div>';
            }
        },
        load: function (query, callback) {
            // if (!query.length) return callback();
            $.ajax({
                url: $('#gene-select').attr('data-url-genes'),
                dataType: 'json',
                error: function () {
                    callback();
                },
                success: function (res) {
                    callback(res);
                    if (localStorage.hasOwnProperty("gene")) {
                        var gene = localStorage['gene']
                        localStorage.removeItem('gene');
                        $gene_select[0].selectize.setValue(gene);
                        first_load = false;
                    }  else if (first_load) {
                        $gene_select[0].selectize.setValue(res[0]['gene_symbol']);
                        first_load = false;
                    } 
                }
            });
        },
        persist: true,
    });
    
        
    
      

    // CHANGE LINKS FOR APPYTERS/ARCHS4/GTEx DYNAMICALLY   
    
    $('#appyter-home').click( async function() {  
        clear_home()
        if ($("#gene-select").val()) {
            var inputvalue = $("#gene-select").val();
            var isChecked= document.getElementById("species-val").checked;
            if (!isChecked) {
                var species = 'Human';
                var arg = 'human_gene';
            } else {
                var species = 'Mouse';
                var arg = 'mouse_gene';
            }
            const formData = new FormData()
            formData.append('species_input', species)
            formData.append(arg, inputvalue)

            document.getElementById("volcano-loading").innerHTML = "<div class='loader mb-2' style='left: 48%; position: relative;'></div>";
            
            var res = await fetch("https://appyters.maayanlab.cloud/Gene_Expression_T2D_Signatures/", {
                method: "POST",
                headers: {
                    'Accept': 'application/json',
                },
                body: formData,
            })

            const id = await res.json()

            const final_url = "https://appyters.maayanlab.cloud/Gene_Expression_T2D_Signatures/" + id.session_id

            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-4 mb-1' onclick='clear_home();'> Clear Results </button> </a>"
            appyter_button = `<a id="appyter-home" href="${final_url}" target='_blank'><button type="button"
                                class="btn btn-primary btn-group-sm mt-3 mb-2">
                                <span id="appyter-action" class="ml-3">Open in</span>
                                <img src="/static/img/appyters_logo.svg" class="img-fluid mr-3" style="width: 120px" alt="Appyters">
                                </button>
                                </a>`
            var gene = inputvalue;
            var jsonData = {};
            species = species.toLowerCase();
            jsonData["gene"] = gene;
            jsonData["species"] = species;
            $.ajax({
                url: "api/volcano",
                type: "POST",
                dataType: 'json',          
                data: jsonData,
                success: function(jdata) {
                    var plot = jdata['plot']
                    var tables = jdata['tables']
                    var micro = jdata['micro']
                    //document.getElementById("buttons").innerHTML += `<div class='row'><div class= 'col text-right'>${appyter_button}</div><div class='col text-left'>${clear_button}</div></div>`
                    document.getElementById("buttons").innerHTML += `<div class='row text-center justify-content-center'>${clear_button}</div>`
                    window.Bokeh.embed.embed_item(plot)
                    document.getElementById("volcano-loading").innerHTML = "";
                    
                    
                        
                    var dir = "up";
                    var titleRNA = `Top ${species} RNA-seq signatures where ${inputvalue} is ${dir}-regulated`
                    var titlemicro = `Top ${species} microarray signatures where ${inputvalue} is ${dir}-regulated`
                    gen_table('', 'faketable', '', false)
                    gen_table(tables[0], `${species}_up`, titleRNA, gene)
                    if (micro) gen_table(tables[2], `${species}_micro_up`, titlemicro, gene)
                    dir = "down";
                    var titleRNA = `Top ${species} RNA-seq signatures where ${inputvalue} is ${dir}-regulated`
                    var titlemicro = `Top ${species} microarray signatures where ${inputvalue} is ${dir}-regulated`
                    gen_table(tables[1], `${species}_down`, titleRNA, gene)
                    if (micro) gen_table(tables[3], `${species}_micro_down`, titlemicro, gene)


                }
            });
            
        } else {
           alert("Please select valid gene symbol")
        }
    });

    $('#appyter-url1').click( async function() {  
        var selectize = $(`#search1`)[0].selectize;
        var gene = selectize.getValue();
        if (gene) {
            var check_list = await human_list
            
            if (check_list.includes(gene)) {
                var species = 'Human';
                var arg = 'human_gene';
            } else{
                var species = 'Mouse';
                var arg = 'mouse_gene';
            }
            const formData = new FormData()
            formData.append('species_input', species)
            formData.append(arg, gene)
            var res = await fetch("https://appyters.maayanlab.cloud/Gene_Expression_by_Tissue/", {
                method: "POST",
                headers: {
                    'Accept': 'application/json',
                },
                body: formData,
            })

            const id = await res.json()
            window.open("https://appyters.maayanlab.cloud/Gene_Expression_by_Tissue/" + id.session_id, target='_blank')


        } else {
            window.open("https://appyters.maayanlab.cloud/Gene_Expression_by_Tissue/", target='_blank')
        }
    });


    $('#archs-url').click(function() {  

        if ($("#search1").val()) {
            var inputvalue = $("#search1").val();
            
            $('#archs-url').prop('href', "https://maayanlab.cloud/archs4/gene/" +inputvalue + "#tissueexpression");
        } else {
            $('#archs-url').prop('href', "https://maayanlab.cloud/archs4/")
        }
    });

    $('#gtex-url').click(function() {  
        if ($("#search1").val()) {
            var inputvalue = $("#search1").val();
            $('#gtex-url').prop('href', "https://gtexportal.org/home/gene/" +inputvalue + "#geneExpression");
        } else {
            $('#gtex-url').prop('href', "https://gtexportal.org/home/")
        }
    });

    $('#appyter-url2').click(async function() {  
        if ($("#search2").val()) {
            var inputvalue = $("#search2").val();
            var check_list = await human_list
            
            if (check_list.includes(inputvalue)) {
                var species = 'Human';
                var arg = 'human_gene';
            } else{
                var species = 'Mouse';
                var arg = 'mouse_gene';
            }
            const formData = new FormData()
            formData.append('species_input', species)
            formData.append(arg, inputvalue)
            var res = await fetch("https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/", {
                method: "POST",
                headers: {
                    'Accept': 'application/json',
                },
                body: formData,
            })
            const id = await res.json()
            window.open("https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/" + id.session_id, target='_blank')
        } else {
            window.open("https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/", target='_blank');
        }
    });

    $('#appyter-url3').click(async function() {  
        if ($("#search2").val()) {
            var inputvalue = $("#search2").val();
            var check_list = await human_list 
            if (!check_list.includes(inputvalue)) {
                alert('This Appyter only accepts Human gene symbols.')
                return;
            }
            window.open("https://appyters.maayanlab.cloud/L1000_RNAseq_Gene_Search/#/?args.gene=" +inputvalue + "&submit", target='_blank');
        } else {
            window.open("https://appyters.maayanlab.cloud/L1000_RNAseq_Gene_Search/", target='_blank');
        }
    });

    $('#appyter-url4').click(function() {  
        if ($("#search3").val()) {
            var inputvalue = $("#search3").val();
            $('#appyter-url4').prop('href', "https://appyters.maayanlab.cloud/ChEA3_Appyter/#/?args.paste_gene_input=" +inputvalue + "&submit");
        } else {
            $('#appyter-url4').prop('href', "https://appyters.maayanhttps://appyters.maayanlab.cloud/ChEA3_Appyter/")
        }
    });

    // QUERY DIABETES PERTURBATIONS ENRICHR LIBRARY

    $('#diabetesEnrichr-query').click(async function() {  
        var inputvalue = document.getElementById("text-area1").value;
        var desc = document.getElementById("desc1").value;
        var file = document.getElementById("gene-file1").value;
        const section = "gene-file1";

        if (!inputvalue && !file) {
            alert("No genes entered")
            return;
        }

        if (file) {
            inputvalue = await loadFileAsText(section, "\n");
        }


        if (!inputvalue) {
            alert("Check file format!")
            return;
        }

        if (!inputvalue) {
            $("#enrich-res").html("");
            return;
        }



        document.getElementById("enrich-res").innerHTML = "<div class='loader' style='left: 48%; position: relative;'></div>"

        $.ajax({
            url: "getdiabetesenrich",
            type: "POST",
            data: {genelist:inputvalue, description: desc}
        }).done(function(response) {

            const data = response['data']['Diabetes_Perturbations_GEO_2022'];



            if (!data) {
                $("#enrich-res").html("<p class='text-center'> No data found </p>");
                return;
            }

            
            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#enrich-res\").html(\"\");'> Clear Results </button> </a>"


            var tabletext = "<table id='table-enrichr' class='styled-table' style:'width=100%; vertical-align:top;'><thead><tr><th>Rank</th><th>Term name</th><th>P-value</th><th>Z-score</th><th>Combined score</th><th>Overlapping genes</th><th>Adjusted p-value</th></tr><tbody>";

            for (var k = 0; k < data.length; k++) {
                tabletext += "<tr><td>" + data[k][0] + "</td><td>"+ data[k][1] +"</td><td>" + Number(data[k][2]).toPrecision(4) + "</td><td>"+Number(data[k][3]).toPrecision(4) +"</td><td>" + Number(data[k][4]).toPrecision(4) + "</td><td>"
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


            $(document).ready(function(){
                $('#table-enrichr').DataTable({
                    dom: 'Bfrtip',
                    buttons: [
                    'copy', {extend: 'csv', title: `${desc}-Diabetes-Perturbations-Enrichr-res`}
                ]});
                    
            });

            document.getElementById("enrich-res").innerHTML = tabletext + clear_button;
        });
    });

    
    // QUERY GWAS AND PRODUCE TABLE

    $('#gwas-query').click(function() {  
        var inputvalue = $("#search4").val();
        if (!inputvalue) {
            $("#gwas-res").html("");
            return;
        }

        document.getElementById("gwas-res").innerHTML = "<div class='loader' style='left: 48%; position: relative;'></div>"

        $.ajax({
            url: "getgwas",
            type: "POST",
            data: {gene:inputvalue}
        }).done(function(response) {

            const data = response['GWAS_Catalog'];


            if (data.length === 0) {
                $("#gwas-res").html("<p class='text-center'> No data found </p>");
                return;
            }

            
            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#gwas-res\").html(\"\");'> Clear Results </button> </a>"


            var tabletext = "<table id='table-gwas' class='styled-table'><thead><tr><th></th><th>Gene</th><th>Trait</th><th>Count</th></tr><tbody>";
            for (var k = 0; k < data.length; k++) {
                tabletext += "<tr><td>"+(k+1)+"</td><td><a href='https://www.ebi.ac.uk/gwas/genes/"+ data[k]['gene']+ "' target='_blank'>"+data[k]['gene']+"</a></td><td><a href='"+ data[k]['mapped_trait_link']+ "' target='_blank'>"+data[k]['trait']+"</a></td><td>"+data[k]['count']+"</td></tr>";
            }
            tabletext += "</tbody></table>";


            $(document).ready(function(){
                $('#table-gwas').DataTable({
                    dom: 'Bfrtip',
                    buttons: [
                    'copy', {extend: 'csv', title: `${inputvalue}-gwas-res`}
                ]});
            });

            document.getElementById("gwas-res").innerHTML = tabletext + clear_button;
        });
    });


    // QUERY ARCHS4 AND CREATE TABLE

    function loadCorrelation(gene){
        $("archs4-res").html("");

        document.getElementById("archs4-res").innerHTML = "<div class='loader' style='left: 48%; position: relative;'></div>"
        var jsonData = {};

        jsonData["id"] = gene;
        jsonData["count"] = 101;
        $.ajax({
            type: "POST",
            url: "https://maayanlab.cloud/matrixapi/coltop",
            contentType: "application/json; charset=utf-8",
            dataType: "json",
            data: JSON.stringify(jsonData),
            success: function(jdata) {
                var data = jdata;

                var genesym = data["rowids"];
                var correlation = data["values"];

                if (!(genesym)) {
                    $("#archs4-res").html("<p class='text-center'> No data found </p>");
                    return;
                }
                const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#archs4-res\").html(\"\");'> Clear Results </button> </a>"
                //3) separate them back out:
                genes = "";
                var tabletext = "<table id='tablecor' class='styled-table'><thead><tr><th>Rank</th><th>Gene Symbol</th><th>Pearson Correlation</th></tr><tbody>";
                for (var k = 1; k < genesym.length; k++) {
                    tabletext += "<tr><td>"+k+"</td><td><a href=\"https://maayanlab.cloud/archs4/gene/"+genesym[k]+"\" target=\"_blank\">"+genesym[k]+"</a></td><td>"+Number(correlation[k]).toPrecision(4)+"</td></tr>";
                    genes = genes+genesym[k]+"\n";
                }
                tabletext += "</tbody></table>";

                $(document).ready(function(){
                    $('#tablecor').DataTable({
                        dom: 'Bfrtip',
                        buttons: [
                        'copy', {extend: 'csv', title: `${gene}-archs4-corr`}
                    ]});
                });

                document.getElementById("archs4-res").innerHTML = tabletext + clear_button;

            },
            error: function (xhr, textStatus, errorThrown) {
            }
        });
    }

    $('#archs4-query').click(function() {  
        var inputvalue = $("#search5").val();
        if (!inputvalue) {
            $("#archs4-res").html("");
            return;
        }
        loadCorrelation(inputvalue);
    });

    // QUERY KOMP AND CREATE TABLE

    $('#komp-query').click(function() {  
        var inputvalue = $("#search6").val();
        if (!inputvalue) {
            $("#komp-res").html("");
            return;
        }
        document.getElementById("komp-res").innerHTML = "<div class='loader' style='left: 48%; position: relative;'></div>"

        $.ajax({
            url: "getkomp",
            type: "POST",
            data: {gene:inputvalue}
        }).done(function(response) {

            const data = response['data'];

            if (data.length === 0) {
                $("#komp-res").html("<p class='text-center'> No data found </p>");
                return;
            }
            //const toggle = "<a id='toggle-vis' data-column='3'> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='hideEvidence();'> Toggle Evidence Column </button> </a>"
            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#komp-res\").html(\"\");'> Clear Results </button> </a>"
            var tabletext = "<table id='table-pheno' class='styled-table'><thead><tr><th>Gene</th><th>Phenotype</th><th>PM ID</th><th>Comments</th></tr><tbody>";
            for (var k = 0; k < data.length; k++) {
                tabletext += "<tr><td><a href='http://www.informatics.jax.org/marker/" + data[k]['OntologyAnnotation.subject.primaryIdentifier'] + "' target='_blank'>"+data[k]['OntologyAnnotation.subject.symbol']+"</a></td>"
                tabletext += "<td><a href='http://www.informatics.jax.org/vocab/mp_ontology/" + data[k]['OntologyAnnotation.ontologyTerm.identifier'] + "' target='_blank'>" + data[k]['OntologyAnnotation.ontologyTerm.name']+ "</a></td>"
                tabletext += "<td><a href='https://pubmed.ncbi.nlm.nih.gov/" + data[k]['OntologyAnnotation.evidence.publications.pubMedId'] + "' target='_blank'>" + data[k]['OntologyAnnotation.evidence.publications.pubMedId']+ "</a></td>"
                tabletext += "<td>" + data[k]['OntologyAnnotation.evidence.comments.description'] +"</td></tr>"
            }
            tabletext += "</tbody></table>";

            $(document).ready(function(){
                $('#table-pheno').DataTable({
                    dom: 'Bfrtip',
                    buttons: [
                    'copy', {extend: 'csv', title: `${inputvalue}-mgi-pheno`}
                ]});
            });

            
            document.getElementById("komp-res").innerHTML = tabletext + clear_button;

        });
    });


    // QUERY ENRICHR AND FORMAT RESULTS

    $('#enrichr-query').click(function() {  
        var inputvalue = $("#search3").val();
        if (!inputvalue) {
            document.getElementById("tf-res").innerHTML = "<p class='text-center'> No data found </p>";
            return;
        }
        document.getElementById("tf-res").innerHTML = "<div class='loader' style='left: 48%; position: relative;'></div>"

        $.ajax({
            url: "gettfs",
            type: "POST",
            data: {gene:inputvalue}
        }).done(function(response) {


            const data = response['data'];

            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#tf-res\").html(\"\");'> Clear Results </button> </a>"

            if (data.length === 0) {
                document.getElementById("tf-res").innerHTML = "<p class='text-center'> No data found </p>";
                return;
            }

            var selecter = `<div class='text-center'><p>Select from one the annotated libraries: </p><select class="m-2 libpicker" data-style="btn-primary" onchange="on_change(this)" data-width="500px">`
            for (var i = 0; i < data.length; i++) {
                var lib = data[i]['name'];
                var libDisplay = lib.replaceAll("_"," ")
                selecter += `<option value=${lib}>${libDisplay}</option>`;
            }
            selecter += `</select></div>`

            

            for (var i = 0; i < data.length; i++) {
                var lib = data[i]['name'];
                var sentence = data[i]['format'];

                var res_html = `<div id="${lib}" style="display:none;"><table id='table-enrichr' class='styled-table table-enrichr'><thead><tr><th></th></tr><tbody>`
                

                for (j = 0; j < data[i]['tfs'].length; j++) {
                    var tf = data[i]['tfs'][j]
             
                    tf_sentence = sentence.replace('{0}', inputvalue).replace('{1}', tf)
                    
                    res_html += `<tr><td> ${tf_sentence} </td></tr>`

                }
                
                res_html += `</tbody></table></div>`

                selecter += res_html

            }
            
            selecter += `<script>
                        
                        </script>`

            
            $(document).ready(function(){

                $(`.table-enrichr`).DataTable({
                    dom: 'Bfrtip',
                    buttons: [
                    'copy', {extend: 'csv', title: `${inputvalue}-enrichr-tfs`}
                ]});
            });
            
        
            document.getElementById("tf-res").innerHTML = selecter + clear_button;
            document.getElementById(data[0]['name']).style.display = 'block';

        });
    });

    
    


    // Configure dropdown menu

    $('.dropdown-menu a.dropdown-toggle').on('click', function(e) {
        if (!$(this).next().hasClass('show')) {
          $(this).parents('.dropdown-menu').first().find('.show').removeClass('show');
        }
        var $subMenu = $(this).next('.dropdown-menu');
        $subMenu.toggleClass('show');
      
      
        $(this).parents('li.nav-item.dropdown.show').on('hidden.bs.dropdown', function(e) {
          $('.dropdown-submenu .show').removeClass('show');
        });
      
      
        return false;
    });


    //////////////////////////////////
    /// Boxplot
    //////////////////////////////////

    // 1. Plotting function

    function boxplot() {

        // Change Status
        $('#boxplotloader').addClass('loader');

        // Gene
        var gene_symbol = $('#gene-select').val();
        
        // Conditions FOR SINGLE CELL NEED TO CHANGE BUT FOR CLUSTERS NOW WHAT ARE THE SELECTED CLUSTERS
        var conditions = [];
        $('.condition-btn.plotted').each(function() { conditions.push($(this).attr('data-group_label')) }); conditions

        // AJAX Query
        $.ajax({
            url: $('#boxplot').attr('data-url-plot'), //"{{ url_for('plot_api') }} + "/" + $('#boxplot').attr('data-geo-acc'),
            method: 'post',
            data: JSON.stringify({'gene_symbol': gene_symbol, 'conditions': conditions}),
            contentType: 'application/json',
            dataType: 'json',
            error: function () {
                "<p> Gene not found <\p>"
            },
            success: function (res) {
                $('#boxplotloader').removeClass('loader');
                layout= {
                    plot_bgcolor: "#00FFFFFF",
                    paper_bgcolor:"#00FFFFFF"
                }
                Plotly.newPlot('boxplot', res['data'], res['layout'], config={responsive: true}); // maybe plotly.react will be faster here
            }
        });

    };


    // 2. Listeners
    // Gene
    var is_gse = currURL.filter(x => x.includes('GSE'))
    if (is_gse.length > 0) {
        var boxplot_selectize = $gene_select[0].selectize;
        boxplot_selectize.on('change', function(value) {
        boxplot();
        })
    }
    
   
    /* $('#generate-plot').on('click', function(evt) {
        boxplot();
    }) */
    
    // Conditions
    $('.condition-btn').on('click', function(evt) {
        $(this).toggleClass('plotted'); // making a specific button plotted or not
        boxplot();
    })

    // 3. Plot
    // boxplot(); // for initial plotting?


    // ADD LIST TO ENRICHR and Redirect

    $('#enrichr').click(async function() {  
        var inputvalue = document.getElementById("text-area1").value;
        var desc = document.getElementById("desc1").value;
        var file = document.getElementById("gene-file1").value;
        const section = "gene-file1";

        if (!inputvalue && !file) {
            alert("No genes entered")
            return;
        }

        if (file) {
            inputvalue = await loadFileAsText(section, "\n");
        }


        if (!inputvalue) {
            alert("Check file format!")
            return;
        }


        if (!desc) {
            enrich({list: inputvalue, popup:true});
        } else {
            enrich({list: inputvalue, description: desc, popup:true});
        }

    });

    // TAKE FILE OR TEXT INPUT AND OPEN IN KEA3

    $('#kea3').click(async function() {  
        var inputvalue = document.getElementById("text-area2").value;
        var file = document.getElementById("gene-file2").value;
        var section = "gene-file2";

        if (inputvalue) {
            inputvalue = inputvalue.split("\n").join("%0A")
        }

        if (file) {
            inputvalue = await loadFileAsText(section, "%0A");
        }

        if (!inputvalue) {
            $('#kea3').prop('href', "https://appyters.maayanlab.cloud/KEA3_Appyter/")
            return;
        }

        $('#kea3').prop('href', "https://appyters.maayanlab.cloud/KEA3_Appyter/#/?args.Input%20gene/protein%20list=" +inputvalue + "&submit");

    });

    // TAKE FILE OR TEXT INPUT AND OPEN IN CHEA3

    $('#chea3').click(async function() {  
        var inputvalue = document.getElementById("text-area2").value;
        var file = document.getElementById("gene-file2").value;
        var section = "gene-file2";

        if (inputvalue) {
            inputvalue = inputvalue.split("\n").join("%0A")
        }

        if (file) {
            inputvalue = await loadFileAsText(section, "%0A");
        }

        if (!inputvalue) {
            $('#chea3').prop('href', "https://appyters.maayanlab.cloud/ChEA3_Appyter/")
            return;
        }

        $('#chea3').prop('href', "https://appyters.maayanlab.cloud/ChEA3_Appyter/#/?args.paste_gene_input=" +inputvalue + "&submit");

    });

    $('#examplefill1').click(function() {
        fillSet('text-area1', 'desc1', 1)
    });

    $('#examplefill2').click(function() {
        fillSet('text-area2', '', 2)
    });

    $('#examplefill3').click(function() {
        fillSet('text-area3', '', 3)
    });


    // PRODUCE INPUT HTML FOR SINGLE GENE SET
    
    function getSingleEntry(num, genecount) {
        var single_entry = `<div class="col-6 col-md-10 col-sm-7 col-lg-10 text-center mr-3">
                  <textarea name="list" rows="8" id="text-area${num}"
                    placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                    onkeyup="geneCount($(this).val(), 3)" onchange="geneCount($(this).val(), 3)"
                    onfocus="geneCount($(this).val(), 3)"></textarea>
                  <div class="mt-1">
                    <span id="gene-count${genecount}"> 0 </span> gene(s) entered
                  </div>
                  <div class="text-center">
                    <a id="examplefill3" style="color: rgb(10, 13, 149)">Try an example gene set</a>
                  </div>
                </div>
                <div class="flex-column justify-content-center text-center m-2" style="overflow-y: visible !important;">
                  <p>
                    File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line
                  </p>
                  <form action="/action_page.php">
                    <input type="file" id="gene-file${num}" name="filename">
                  </form>

                </div>`;
        return single_entry;
    }

    // PRODUCE INPUT HTML FOR MULTIPLE GENE SETS

    function getMultipleEntries(num, genecount) {
        var multiple_entries = `<div class="col-6 col-md-10 col-sm-7 col-lg-10 text-center mr-3">
            <textarea name="list" rows="8" id="text-area${num}-up" placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box" onkeyup="geneCount($(this).val(), ${genecount})" onchange="geneCount($(this).val(),  ${genecount})" onfocus="geneCount($(this).val(),  ${genecount})"></textarea>
            <div class="mt-1">
                <span id="gene-count${genecount}"> 0 </span> UP gene(s) entered
            </div>
            <div class="text-center">
                <a id="examplefill3-up" onclick="fillSet('text-area3-up', '', 3)" style="color: rgb(10, 13, 149)">Try an example gene set</a>
            </div>
            </div>
            <div class="flex-column justify-content-center text-center m-2" style="overflow-y: visible !important;">
            <p>
                File formats excepted: csv, tsv, txt file with Entrez gene symbols on each line
            </p>
            <form action="/action_page.php">
                <input type="file" id="gene-file${num}-up" name="filename">
            </form>
            
            </div>
            <div class="col-6 col-md-10 col-sm-7 col-lg-10 text-center mr-3">
            <textarea name="list" rows="8" id="text-area${num}-down" placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box" onkeyup="geneCount($(this).val(), ${genecount + 1})" onchange="geneCount($(this).val(), ${genecount + 1})" onfocus="geneCount($(this).val(), ${genecount + 1})"></textarea>
            <div class="mt-1">
                <span id="gene-count${genecount + 1}"> 0 </span> DOWN gene(s) entered
            </div>
            <div class="text-center">
                <a id="examplefill3-down" onclick="fillSet('text-area3-down', '', 4)" style="color: rgb(10, 13, 149)">Try an example gene set</a>
            </div>
            </div>
            <div class="flex-column justify-content-center text-center m-2" style="overflow-y: visible !important;">
            <p>
                File formats excepted: csv, tsv, txt file with Entrez gene symbols on each line
            </p>
            <form action="/action_page.php">
                <input type="file" id="gene-file${num}-down" name="filename">
            </form>
            
            </div>`;
        return multiple_entries;
    }

    // SWITCH INPUT FOR SIGCOM LINCS to UP DOWN OR SINGLE GENE SET

    $('#sigcom_entries_button').click( function() {  
        
        var button = document.getElementById("sigcom_entries_button");
        var mode = button.value


        if (mode === "single"){
            document.getElementById("sigcom_entries").innerHTML = getMultipleEntries(3, 3);
            button.innerText = "Use Single Gene Set"
            button.value = "double"
        } else {
            document.getElementById("sigcom_entries").innerHTML = getSingleEntry(3, 3);
            button.innerText = "Use Up/Down Gene Sets"
            button.value = "single"
        }
        
    });

    // READ IN TEXT/FILE INPUT AND OPEN IN SIGCOM LINCS

    $('#sigcoms').click(async function() {  
        var button = document.getElementById("sigcom_entries_button");


        if (button.textContent === "Use Up/Down Gene Sets") {

            var inputvalue = document.getElementById("text-area3").value;
            var file = document.getElementById("gene-file3").value;
            var section = "gene-file3";

            if (inputvalue) {
                inputvalue = inputvalue.split("\n")
            }

            if (file) {
                inputvalue = await loadFileAsText(section, "\t").split("\t");
            }

            if (inputvalue.length == 0 ) {
                alert('Please enter a gene set')
            }

            var genes = JSON.stringify({'genes': [inputvalue]});

            $.ajax({
                url: "getsigcom",
                contentType: 'application/json',
                type: "POST",
                dataType: 'json',
                data: genes
            }).done(function(response) {
    
                const url = response['url'];
                window.open(url, '_blank');

            });

        } else {
            var inputvalueUp = document.getElementById("text-area3-up").value;
            var inputvalueDown = document.getElementById("text-area3-down").value;

            var fileUp = document.getElementById("gene-file3-up").value;
            var sectionUp = "gene-file3";
            var fileDown = document.getElementById("gene-file3-down").value;
            var sectionDown = "gene-file3";

            if (inputvalueUp) {
                inputvalueUp = inputvalueUp.split("\n")
            }
            if (inputvalueDown) {
                inputvalueDown = inputvalueDown.split("\n")
            }

            if (fileUp) {
                inputvalueUp = await loadFileAsText(section, "\t").split("\t");
            }

            if (fileDown) {
                inputvalueDown = await loadFileAsText(section, "\t").split("\t");
            }

            if (inputvalueUp.length == 0 || inputvalueDown.length == 0) {
                alert('Please enter an up and down gene set')
            }

            var genes = JSON.stringify({'genes': [inputvalueUp, inputvalueDown]});

            $.ajax({
                url: "getsigcom",
                contentType: 'application/json',
                type: "POST",
                dataType: 'json',
                data: genes
            }).done(function(response) {
    
                const url = response['url'];
                window.open(url, '_blank');
            });

        }
        
        
    });

    // MAKE STUIDIES TABLE A DataTable

    $('#studies-table').DataTable({
        dom: 'Bfrtip',
        buttons: [
        'copy', {extend: 'csv', title: `D2H2-studies-table`}
        ],
        responsive: true
    });


    // OPEN CUSTOMIZED WORKFLOW DEPENDING ON SELECTION IN SCG

    $('#scg-link').click( function() {

        var workflow = document.getElementById("workflow").value;
        

        if (workflow) {

            $('#scg-link').prop('href', workflow);
            return;
        }
        $('#scg-link').prop('href', "https://scg.maayanlab.cloud/");
    });

    // OPEN TO DEG IN BULK RNA SEQ ANALYSIS

    $('#dgea-button').on('click', async function() {
        var control_condition = $('#condition-select-control').val();
        var perturb_condition = $('#condition-select-perturb').val();
        var species = document.getElementById("species").innerText

        if (!control_condition || !perturb_condition ) {
            alert("Please select both conditions")
            return;
        }

        if (control_condition == perturb_condition) {
            alert("Please select two different conditions")
            return;
        }

        document.getElementById("dgea-loading").innerHTML = "<div class='loader justify-content-center'></div>";

        var gse = document.getElementById("gse").innerText

        var gsedata = JSON.stringify({'gse': gse, 'control': control_condition, 'perturb': perturb_condition, 'species': species});


        
        $.ajax({
            url: "api/data",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done(async function(response) {

            meta = response['meta']
            expression = response['expression']

            meta_data_file = gse + '_Metadata.txt'
            rnaseq_data_filename = gse + '_Expression.txt'

            const blob_meta = new Blob([meta]);
            const blob_rna = new Blob([expression]);
            
            
            const formData = new FormData()
            formData.append('control_name', control_condition)
            formData.append('meta_class_column_name', 'Condition')
            formData.append('meta_data_filename', blob_meta, meta_data_file)
            formData.append('rnaseq_data_filename', blob_rna, rnaseq_data_filename)
            
            var res = await fetch("https://appyters.maayanlab.cloud/Bulk_RNA_seq/", {
                method: "POST",
                headers: {
                    'Accept': 'application/json',
                },
                body: formData,
            })


            const id = await res.json()

            document.getElementById("dgea-loading").innerHTML = "";

            var download_link = "https://appyters.maayanlab.cloud/Bulk_RNA_seq/" + id.session_id + "/DEG_results_" + control_condition.split(" ").join('%20') + "%20vs.%20"+ perturb_condition.split(" ").join('%20') + ".csv"


            const final_url = "https://appyters.maayanlab.cloud/Bulk_RNA_seq/" + id.session_id + "/#differential-gene-expression"
            window.open(final_url, '_blank')

        });
    });

    // PERFORM DIFFERENTIAL GENE ANALYSIS AND CREATE RELEVANT TABLE
    $('#dge-button').on('click', async function() {
        clear_dge()
        var control_condition = $('#condition-select-control').val();
        var perturb_condition = $('#condition-select-perturb').val();

        if (!control_condition || !perturb_condition ) {
            alert("Please select both conditions")
            return;
        }

        if (control_condition == perturb_condition) {
            alert("Please select two different conditions")
            return;
        }

        document.getElementById("dge-loading").innerHTML = "<div class='loader justify-content-center'></div> <div class='text-center'><p> It may take 2-3 minutes to compute DEGs</p><div>";

        var gse = document.getElementById("gse").innerText
        var species = document.getElementById("species").innerText
        var method = document.getElementById("method").value
         
        var logCPM = document.getElementById("logCPM").checked
        var log = document.getElementById("log").checked
        var q = document.getElementById("q").checked
        var z = document.getElementById("z").checked
        var norms = {'logCPM': logCPM, 'log': log, 'q': q, 'z': z}


        var gsedata = JSON.stringify({'gse': gse, 'control': control_condition, 'perturb': perturb_condition, 'method': method, 'species': species, 'norms': norms});
        $.ajax({
            url: "dgeapi",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done(async function(response) {
            document.getElementById("dge-loading").innerHTML = "";
            var plot = response['plot']
            var table = response['table']

            var rows = table.split('\n').slice(1, -1);
            clear_dge()
            window.Bokeh.embed.embed_item(plot)

            var table_id = 'dge-table'
            var tabletext = `<table id='${table_id}' class='styled-table'><thead><tr>`
            const clear_button = "<div class='mx-auto justify-content-center text-center'><button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_dge();'> Clear Results </button></div>"
            var url = currURL.filter(x => !x.includes('GSE')).join('/') + '/singlegene'
            var pvals;
            if (method === 'limma') {
                tabletext += "<th></th><th>Adj. P Value</th><th>P Value</th><th>t</th><th>AvgExpr</th><th>logFC</th><th>B</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>`  + vals[0] + "</a></td><td>" + Number(vals[5]).toPrecision(4)+"</td><td>" + Number(vals[4]).toPrecision(4) +"</td><td>"+ Number(vals[3]).toPrecision(4) + "</td><td>"+ Number(vals[2]).toPrecision(4) + "</td><td>"+ Number(vals[1]).toPrecision(4) + "</td><td>"+ Number(vals[6]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()

            } else if (method === 'edgeR') {
                tabletext += "<th></th><th>PValue</th><th>logCPM</th><th>logFC</th><th>FDR</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>`  +vals[0] + "</a></td><td>" + Number(vals[3]).toPrecision(4)+"</td><td>" + Number(vals[2]).toPrecision(4) +"</td><td>"+ Number(vals[1]).toPrecision(4) + "</td><td>"+ Number(vals[4]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext 
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(1).data()

            } else if (method === 'DESeq2') {
                tabletext += "<th></th><th>Adj. P-value</th><th>P-value</th><th>lfcSE</th><th>stat</th><th>baseMean</th><th>log2FC</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` +vals[0] +"</a></td><td>" + Number(vals[6]).toPrecision(4)+"</td><td>" + Number(vals[5]).toPrecision(4) +"</td><td>"+  Number(vals[3]).toPrecision(4)  + "</td><td>"+  Number(vals[4]).toPrecision(4)  + "</td><td>" +  Number(vals[1]).toPrecision(4)  + "</td><td>" +  Number(vals[2]).toPrecision(4)  + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });    
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
            }
            var genes = table.column(0).data()

            
            genes = genes.map(x => x.replace(/<\/?[^>]+(>|$)/g, ""))

            var genelist_buttons = 
            `<div class="row justify-content-center mx-auto text-center">
            <div class="h7">Submit the top</div>
            <input class="" id='numgenes' type='number' step='1' value='100' pattern='[0-9]' min='1' class='m-2' style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;'/>
            <div class="h7">  differentially expressed genes with a 
            <select id='col-to-use' style='margin-left: 10px; margin-right: 10px;'><option value='pval'>p-value</option><option value='adjpval'>adjusted p-value</option></select></div>
            <div class=" h7">  less than  </div>
            <input class='' id='signifigance' type='number' step='.01' value='.05' max='1' style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;"'/>
            <div class="h7">  to</div>
            </div>
            <div class="row justify-content-center mx-auto text-center">
            <button class="btn btn-primary btn-group-sm m-2" onclick="submit_geneset('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}')">Gene Set Queries</button>
            <button class="btn btn-primary btn-group-sm m-2" onclick="submit_geneset_home('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${name}')">Diabetes Gene Set Library</button>
            </div>
            <div class="row justify-content-center mx-auto text-center">
            <button type="button" class="btn btn-primary btn-group-sm m-2" onclick="filter_and_submit_to_enrichr('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${gse}-${control_condition}-vs-${perturb_condition}')"> Enrichr
            <img src="/static/img/enrichrlogo.png" class="img-fluid mr-3" style="width: 45px" alt="Enrichr">
            </button>
            <button type="button" class="btn btn-primary btn-group-sm m-2" onclick="filter_and_submit_to_kg('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${gse}-${control_condition}-vs-${perturb_condition}')"> Enrichr-KG
            <img src="/static/img/enrichr-kg.png" class="img-fluid mr-3" style="width: 45px" alt="Enrichr">
            </button>
            </div>
            `
                
            document.getElementById("geneset-buttons").innerHTML = (clear_button + genelist_buttons)
        })
        .fail(function (jqXHR, textStatus, errorThrown) {
            alert('An internal server error occured, please try again')
            document.getElementById("dge-loading").innerHTML = "";
        })
    });

   

    // PERFORM DIFFERENTIAL GENE ANALYSIS AND CREATE RELEVANT TABLE
    $('#dge-button-single').on('click', async function() {
        clear_dge_single()

        document.getElementById("dge-loading").innerHTML = "<div class='loader justify-content-center'></div>";

        var gse = document.getElementById("singlegse").innerText
        var species = document.getElementById("species").innerText
        var method = document.getElementById("method").value
        var condition_group = document.getElementById("methodsingle").value
        var diffcluster = document.getElementById("differentialcluster").value
         
        var logCPM = document.getElementById("logCPM").checked
        var log = document.getElementById("log").checked
        var q = document.getElementById("q").checked
        var z = document.getElementById("z").checked
        var norms = {'logCPM': logCPM, 'log': log, 'q': q, 'z': z}

        var gsedata = JSON.stringify({'gse': gse, 'species': species, 'conditiongroup':condition_group, 'method': method, 'norms': norms, 'diffcluster':diffcluster});

        $.ajax({
            url: "/dgeapisingle",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done(async function(response) {
            document.getElementById("dge-loading").innerHTML = "";
            var plot = response['plot']
            var table = response['table']
            var desc = response['description']

            var rows = table.split('\n').slice(1, -1);
            clear_dge_single()
            window.Bokeh.embed.embed_item(plot)

            var table_id = 'dge-table'
            var tabletext = `<table id='${table_id}' class='styled-table'><thead><tr>`
            const clear_button = "<div class='mx-auto justify-content-center text-center'><button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_dge_single();'> Clear Results </button></div>"
            var url = currURL.filter(x => !x.includes('GSE')).join('/') + '/singlegene'
            var pvals;
            if (method === 'limma') {
                tabletext += "<th></th><th>Adj. P Value</th><th>P Value</th><th>t</th><th>AvgExpr</th><th>logFC</th><th>B</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>`  + vals[0] + "</a></td><td>" + Number(vals[5]).toPrecision(4)+"</td><td>" + Number(vals[4]).toPrecision(4) +"</td><td>"+ Number(vals[3]).toPrecision(4) + "</td><td>"+ Number(vals[2]).toPrecision(4) + "</td><td>"+ Number(vals[1]).toPrecision(4) + "</td><td>"+ Number(vals[6]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()

            } else if (method === 'edgeR') {
                tabletext += "<th></th><th>PValue</th><th>logCPM</th><th>logFC</th><th>FDR</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>`  +vals[0] + "</a></td><td>" + Number(vals[3]).toPrecision(4)+"</td><td>" + Number(vals[2]).toPrecision(4) +"</td><td>"+ Number(vals[1]).toPrecision(4) + "</td><td>"+ Number(vals[4]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext 
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(1).data()

            } else if (method === 'DESeq2') {
                tabletext += "<th></th><th>Adj. P-value</th><th>P-value</th><th>lfcSE</th><th>stat</th><th>baseMean</th><th>log2FC</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` +vals[0] +"</a></td><td>" + Number(vals[6]).toPrecision(4)+"</td><td>" + Number(vals[5]).toPrecision(4) +"</td><td>"+  Number(vals[3]).toPrecision(4)  + "</td><td>"+  Number(vals[4]).toPrecision(4)  + "</td><td>" +  Number(vals[1]).toPrecision(4)  + "</td><td>" +  Number(vals[2]).toPrecision(4)  + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });    
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
            }else if (method === 'wilcoxon') {
                tabletext += "<th></th><th>Adj. P-value</th><th>P-value</th><th>logfoldchanges</th><th>scores</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>`  +vals[0] + "</a></td><td>" + Number(vals[4]).toPrecision(4)+"</td><td>" + Number(vals[3]).toPrecision(4) +"</td><td>"+ Number(vals[2]).toPrecision(4) + "</td><td>"+ Number(vals[1]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });    
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
            }
            var genes = table.column(0).data()

            
            genes = genes.map(x => x.replace(/<\/?[^>]+(>|$)/g, ""))

                
            document.getElementById("geneset-buttons").innerHTML = (clear_button)
            document.getElementById("enrichment-area").innerHTML  = `
            <div class="h4 pl-2 mt-4 mb-4 text-center">Enrichment Analysis for Highly Expressed Genes in ${desc}</div>
            <div class="row justify-content-center mx-auto text-center">
            <div class="h7">Submit the top</div>
            <input class="" id='numgenes' type='number' step='1' value='100' pattern='[0-9]' min='1' class='m-2' style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;'/>
            <div class="h7">  differentially expressed genes with a
            <select id='col-to-use' style='margin-left: 10px; margin-right: 10px;'><option value='pval'>p-value</option><option value='adjpval'>adjusted p-value</option></select></div>
            <div class=" h7">  less than  </div>
            <input class='' id='signifigance' type='number' step='.01' value='.05' max='1' style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;"'/>
            <div class="h7">  to</div>
            </div>
            <div class="row justify-content-center mx-auto text-center">
            <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3" onclick="filter_and_submit_to_enrichr('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', 'Highly Expressed Genes in ${desc}')"> Submit to Enricher
            <img src="/static/img/enrichrlogo.png" class="img-fluid mr-3" style="width: 45px" alt="Enrichr">
            </button>
            </div>
            `
        })
        .fail(function (jqXHR, textStatus, errorThrown) {
            alert('An internal server error occured, please try again')
            document.getElementById("dge-loading").innerHTML = "";
        })
    });





    //////////////////////////////////
    /// Control Condition Selection
    //////////////////////////////////
    $('#condition-select-control').selectize({
        preload: true,
        valueField: 'Condition',
        labelField: 'Condition',
        searchField: 'Condition',
        render: {
            option: function (item, escape) {
                return '<div class="pt-2 light">'+item.Condition+'</div>';
            }
        },
        load: function (query, callback) {
            // if (!query.length) return callback();
            $.ajax({
                url: $('#condition-select-control').attr('data-url-control'),//"{{ url_for('conditions_api', geo_accession=geo_accession) }}",
                dataType: 'json',
                error: function () {
                    callback();
                },
                success: function (res) {
                    callback(res);
                }
            });
        }
    });

     //////////////////////////////////
    /// Perturb Condition Selection
    //////////////////////////////////
    $('#condition-select-perturb').selectize({
        preload: true,
        valueField: 'Condition',
        labelField: 'Condition',
        searchField: 'Condition',
        render: {
            option: function (item, escape) {
                return '<div class="pt-2 light">'+item.Condition+'</div>';
            }
        },
        load: function (query, callback) {
            // if (!query.length) return callback();
            $.ajax({
                url: $('#condition-select-perturb').attr('data-url-perturb'),//"{{ url_for('conditions_api', geo_accession=geo_accession) }}",
                dataType: 'json',
                error: function () {
                    callback();
                },
                success: function (res) {
                    callback(res);
                }
            });
        }
    });

    $('#species-val').on('change', function() {
        if ($(this).is(':checked')) {
            switchStatus = $(this).is(':checked');
            $gene_select[0].selectize.clearOptions();
            $gene_select[0].selectize.load(function(callback) {
            $.ajax({
                url: 'api/genes/mouse',
                dataType: 'json',
                error: function () {
                    callback();
                },
                success: function (res) {

                    callback(res);
                }
            });
        });         
        }
        else {
            switchStatus = $(this).is(':checked');
            $gene_select[0].selectize.clearOptions();
            $gene_select[0].selectize.load(function(callback) {
            $.ajax({
                url: 'api/genes/human',
                dataType: 'json',
                error: function () {
                    callback();
                },
                success: function (res) {

                    callback(res);
                }
            });
        });         
        }
    })


    //This listener is for when we click on a different group of sample individuals to look at for a study within the single cell viewer on the dropdown select menu. 
    $('#methodsingle').on('change', async function() {
        clear_dge_single()
        document.getElementById("change-loading").innerHTML = "<div class='loader justify-content-center'></div>";
        var gse = document.getElementById("singlegse").innerText
        var species = document.getElementById("species").innerText
        var condition_group = document.getElementById("methodsingle").value
        var gsedata = JSON.stringify({'gse': gse, 'species': species, 'conditiongroup':condition_group});
        generate_single_plots()

        $.ajax({
            url: "/getclusterdata",
            contentType: 'application/json',
            type: "POST",
            data: gsedata,
            dataType: 'json'
        }).done(async function(response) {
    
            const classes = response['classes']
            const metadict = response['metadict']

            document.getElementById("singlecell").innerHTML = ''
            document.getElementById("differentialcluster").innerHTML = ''
            tabletext = ''
            diffselectiontext = ''
            for (var k = 0; k < classes.length; k++) {
                    tabletext += `<tr>
                            <td class="p-0 ml-3 border"><button
                                    style="background-color: #ECEFF1; color: black; width: max-content;"
                                    class="btn m-0 rounded-0 py-0 condition-btn active plotted " data-toggle="button"
                                    autocomplete="off"
                                    data-group_label="${classes[k]}">${classes[k]}</button>
                            </td>
    
                            <td class="border" style="padding: 0px 18px; width: max-content;">
                                <button
                                    style="width: max-content;"
                                    class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-center text-center"
                                    data-toggle="collapse" data-target="#samples-${gse}" aria-expanded="false"
                                    aria-controls="samples-${gse}">
                                    <div class="text">${metadict[classes[k]]} Cells
                                    </div>
                                </button>
    
    
                            </td>
                            </tr>`
                    diffselectiontext += `<option value="${classes[k]}">${classes[k]}</option>`
                }
            document.getElementById("singlecell").innerHTML = tabletext
            document.getElementById("differentialcluster").innerHTML = diffselectiontext
            

            //Add the listener for the condition buttons back in this function 
            $('.condition-btn').on('click', function(evt) {
                $(this).toggleClass('plotted'); // making a specific button plotted or not
                boxplot();
            })
            //Fill up the gene selectize with the proper gene list based off each study. 
            $gene_select[0].selectize.clearOptions();
            $gene_select[0].selectize.load(function(callback) {
            $.ajax({
                url: `/api/singlegenes/${gse}/${condition_group}`,
                dataType: 'json',
                error: function () {
                    callback();
                },
                success: function (res) {

                    callback(res);
                    $gene_select[0].selectize.setValue(res[0]['gene_symbol']);
                    $("#boxplot").attr("data-url-plot", `/api/plot_single/${gse}/${condition_group}`)
                    boxplot() 
                }
            });
        });
    
    
        });

    document.getElementById("change-loading").innerHTML = ''
    })
})









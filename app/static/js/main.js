
function createTweetsTable() {
    document.getElementById("tweets-res").innerHTML = "<div class='loader justify-content-center'></div>";
    $.ajax({
        url: "/gettweets",
        type: "POST",
        data: {},
        dataType: 'json',
    }).done(function(response) {

        const data = response['tweets']

        var headers = data[0]

        var tabletext = "<table id='table-twitter' class='styled-table'><thead><tr>"

        headers.forEach(function(header) {
            tabletext += "<th>" + header + "</th>"
        });

        tabletext += "</tr><tbody>"
        
        for (var k = 1; k < data.length; k++) {


            tabletext += "<tr><td>"+ data[k][0]+"</td><td>"+ data[k][1]+ "</td><td>" + data[k][2] + "</td><td>"+ data[k][3] + "</td>"
            tabletext += "<td><a href='" + data[k][4] +"' target='_blank'>"+ 'link' +"</a></td>"

            var analyze = data[k][5].replaceAll('[', '').replaceAll(']', '').replaceAll(' ', '').replaceAll("'","").split(',')
            tabletext += "<td><a href='" + analyze[0] +"' target='_blank'>"+ "<img class='mr-2' src='static/img/d2h2logo.png' style='width: 25px;'/>" +"</a>"
            tabletext += "<a href='" + analyze[1] +"' target='_blank'>"+ "<img class='mr-2' src='static/img/enrichrlogo.png' style='width: 25px;'/>" +"</a>"
            tabletext += "<a href='" + analyze[2] +"' target='_blank'>"+ "<img class='mr-2' src='static/img/harmonizomelogo.png' style='width: 25px;'/>" +"</a>"
            tabletext += "</td></tr>"
            
        }
        
        tabletext += "</tbody></table>";


        $(document).ready(function(){
            document.getElementById("tweets-res").innerHTML = tabletext;
            table = $('#table-twitter').DataTable({"pageLength": 5});
        });

    });
}
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

function submit_geneset(genelist, sigs) {
    genelist = genelist.split(',')
    sigs = sigs.split(',').map(function(item) {
        return parseInt(item, 10);
    });
    var numgenes = document.getElementById('numgenes').value
    var signifigance = document.getElementById('signifigance').value
    var dir = document.getElementById('dir').value
    var genes_valid = []

    for (i=0; i < genelist.length; i++ ){
        if (sigs[i] <= signifigance) {
            genes_valid.push(genelist[i])
        }
    }
    if (dir === 'top') {
        var genes = genes_valid.splice(0, numgenes).join('&')
        
    } else {
        var genes = genes_valid.slice(-numgenes).join('&')
    }


    const currURL = window.location.href.split('/')
    var url = currURL.splice(0, 3).join('/') + '/geneset/' + genes
    console.log(url)
    window.open(url, '_blank')
}

function submit_geneset_home(genelist, sigs, descset) {
    genelist = genelist.split(',')
    sigs = sigs.split(',').map(function(item) {
        return parseInt(item, 10);
    });
    var numgenes = document.getElementById('numgenes').value
    var signifigance = document.getElementById('signifigance').value
    var dir = document.getElementById('dir').value

    for (i=0; i < genelist.length; i++ ){
        if (sigs[i] <= signifigance) {
            genes_valid.append(genes[i])
        }
    }
    if (dir === 'top') {
        var genes = genes_valid.splice(0, numgenes).join('&')
    } else {
        var genes = genes_valid.slice(-numgenes).join('&')
    }
    localStorage.setItem('genes', genes)
    localStorage.setItem('descset', `${descset}-${dir}-${numgenes}`)
    var home = window.location.href.split('/').filter(x => !x.includes('GSE')).join('/')
    window.open(home, '_blank')
}

function clear_dge() {
    document.getElementById("dge-table-area").innerHTML = ""
    document.getElementById("dge-plot").innerHTML = ""
    document.getElementById("dge-loading").innerHTML = ""
    document.getElementById("geneset-buttons").innerHTML = ""
    
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

function loadFileAsText(section, delim){
    
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
        console.log(geneset)
        document.getElementById(this.id).value = geneset;
        geneCount(geneset, this.id[this.id.length -1])
    })
}


function fillSet(id, descid, count_id) {
    $.ajax({
        url: "/getexample",
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

    var currURL = window.location.href.split("/");

    function createResourcesTable() {
        document.getElementById("resources").innerHTML = "<div class='loader justify-content-center'></div>";
        $.ajax({
            url: "/getresources",
            type: "POST",
            data: {},
            dataType: 'json',
        }).done(function(response) {
    
            const data = response['resources']
    
            var headers = data[0]
    
            var tabletext = "<table id='table-resources' class='styled-table'><thead><tr>"
    
            headers.forEach(function(header) {
                if (header != 'URL') {
                    tabletext += "<th>" + header + "</th>"
                }     
            });
            tabletext += "</tr><tbody>"
    
            for (var k = 1; k < data.length; k++) {
                tabletext += "<tr><td><a href='" + data[k][2] +"' target='_blank'>"+data[k][0]+"</a></td><td>"+ data[k][1]+ "</td>"
                if (data[k][3] != 'N/A') {
                    tabletext += "<td><a href='https://pubmed.ncbi.nlm.nih.gov/" + data[k][3] +"' target='_blank'>"+data[k][3]+"</a></td>"
                } else {
                    tabletext += "<td>"+ data[k][3]+ "</td>"
                }
                tabletext += "<td>"+ data[k][4]+ "</td><td>" + data[k][5] + "</td></tr>"
    
            }
    
            tabletext += "</tbody></table>";
    
    
            $(document).ready(function(){
                document.getElementById("resources").innerHTML = tabletext;
                table = $('#table-resources').DataTable({
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: 'D2H2-resourceslist'}
                    ]
                });
            });

        });
    }

    if (currURL[3] == 'resources') {createResourcesTable();}
    if (currURL[3] == 'downloads') {createDownloadsTable();}
    if (currURL[3] == 'scg') {createWorkflowsTable();}

    


    function createWorkflowsTable() {
        document.getElementById("workflows").innerHTML = "<div class='loader justify-content-center'></div>";
        $.ajax({
            url: "/getworkflows",
            type: "POST",
            data: {},
            dataType: 'json',
        }).done(function(response) {
    
            const data = response['workflows']
    
            var headers = data[0]
    
            var tabletext = "<table id='table-workflows' class='styled-table'><thead><tr>"
    
            headers.forEach(function(header) {
                tabletext += "<th>" + header + "</th>"
            });
            tabletext += "</tr><tbody>"
    
            for (var k = 1; k < data.length; k++) {
                tabletext += "<tr><td>"+ data[k][0]+"</td><td>"+ data[k][1]+ "</td>"
                tabletext += `<td><a id="scg-link" href="${data[k][2]}" target="_blank">
                                    <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3">
                                        Open in<img src="/static/img/scglogo.png" class="img-fluid" style="width: 50px" alt="SCG">
                                    </button>
                              </a></td></tr>`
            }
    
            tabletext += "</tbody></table>";
    
    
            $(document).ready(function(){
                document.getElementById("workflows").innerHTML = tabletext;
                table = $('#table-workflows').DataTable();
            });

        });
    }


    function createDownloadsTable() {
        document.getElementById("downloads").innerHTML = "<div class='loader justify-content-center'></div>";
        $.ajax({
            url: "/getdownloads",
            type: "POST",
            data: {},
            dataType: 'json',
        }).done(function(response) {
    
            const data = response['downloads']
    
            var headers = data[0]
    
            var tabletext = "<table id='table-downloads' class='styled-table mb-5'><thead><tr>"
    
            headers.forEach(function(header) {
                tabletext += "<th>" + header + "</th>"
            });
            tabletext += "</tr><tbody>"
    
            for (var k = 1; k < data.length; k++) {
                tabletext += "<tr><td>"+ data[k][0]+"</td><td>"+ data[k][1]+ "</td><td>"+ data[k][2]+ "</td><td>"

                var links = data[k][3].split(',')
                links.forEach(function(link) {
                    var types = link.split('.')
                    var type = types[types.length - 1]
                    if (type === 'f') {
                        type = 'feather';
                    }
                    tabletext += "<a href='" + link + "'><img src='static/img/download.png' alt='' style='width: 12px;'>"+ "<img class='mr-2 ml-1' src='static/img/" + type + ".png' alt='' style='width: 15px;'></a>";
                });
                tabletext += "</td></tr>"
            }
    
            tabletext += "</tbody></table>";
    
    
            $(document).ready(function(){
                document.getElementById("downloads").innerHTML = tabletext;
                table = $('#table-downloads').DataTable();
            });
        });
    }


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
        if (url.split('/')[3].startsWith('GSE')) {
            $("#viewer").addClass('active'); 
            $("#viewer").parents('li').addClass('active');
        }
    });
    

    




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
                        }
                    }
                });
            },
            persist: false,
        });
    
        
    
      

    // CHANGE LINKS FOR APPYTERS/ARCHS4/GTEx DYNAMICALLY   
    
    $('#appyter-home').click( async function() {  
        clear_home()
        if ($("#gene-select").val()) {
            var inputvalue = $("#gene-select").val();
            var isChecked= document.getElementById("species-val").checked;
            if (!isChecked) {
                console.log('human')
                var species = 'Human';
                var arg = 'human_gene';
            } else {
                console.log('mouse')
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
                url: "/api/volcano",
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

    $('#appyter-url1').click(function() {  
        var selectize = $(`#search1`)[0].selectize;
        var gene = selectize.getValue();
        if (gene) {
            $('#appyter-url1').prop('href', "https://appyters.maayanlab.cloud/Gene_Expression_by_Tissue/#/?args.gene=" + gene + "&submit");
        } else {
            $('#appyter-url1').prop('href', "https://appyters.maayanlab.cloud/Gene_Expression_by_Tissue/")
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

    $('#appyter-url2').click(function() {  
        if ($("#search2").val()) {
            var inputvalue = $("#search2").val();
            $('#appyter-url2').prop('href', "https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/#/?args.human_gene=" +inputvalue + "&submit");
        } else {
            $('#appyter-url2').prop('href', "https://appyters.maayanlab.cloud/Gene_Centric_GEO_Reverse_Search/")
        }
    });

    $('#appyter-url3').click(function() {  
        if ($("#search2").val()) {
            var inputvalue = $("#search2").val();
            $('#appyter-url3').prop('href', "https://appyters.maayanlab.cloud/L1000_RNAseq_Gene_Search/#/?args.gene=" +inputvalue + "&submit");
        } else {
            $('#appyter-url3').prop('href', "https://appyters.maayanlab.cloud/L1000_RNAseq_Gene_Search/")
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
            url: "/getdiabetesenrich",
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
                var api = currURL.join('/') + 'singlegene/'
                var api2 = currURL.join('/') + 'geneset/'
                var gene_arr = data[k][5].map(g => `<a href='${api}${g}' target='_blank'>${g}<a/>`);
    
                tabletext += `<button class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-start text-left"
                        data-toggle="collapse" data-target="#genesoverlap-${data[k][0]}" aria-expanded="false"
                        aria-controls="genesoverlap-${data[k][0]}">
                        <div class="text">Show Overlapping Genes</div>
                </button>
                    <div class="collapse" id="genesoverlap-${data[k][0]}">
                        ${gene_arr.join(", ")}
                        <a href="${api2}${data[k][5].join('&')}" target='_blank'>
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
            url: "/getgwas",
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
            url: "/getkomp",
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
            url: "/gettfs",
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
        $('#boxplot').addClass('loading');

        // Gene
        var gene_symbol = $('#gene-select').val();

        if (!gene_symbol) {

            gene_symbol = 'A1CF';

            if(document.getElementById("species")) {
                var species = document.getElementById("species").innerText
            }

            if (species === 'mouse') {
                gene_symbol = '0610007P14Rik';
            }
        }

        // Conditions
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
                $('#boxplot').removeClass('loading');
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
    if (currURL[3].startsWith('GSE')) {
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
    boxplot(); // for initial plotting?



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
                url: "/getsigcom",
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
                url: "/getsigcom",
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

        var gse = currURL[3]
        var gsedata = JSON.stringify({'gse': gse, 'control': control_condition, 'perturb': perturb_condition, 'species': species});


        
        $.ajax({
            url: "/api/data",
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

        document.getElementById("dge-loading").innerHTML = "<div class='loader justify-content-center'></div>";

        var gse = document.getElementById("gse").innerText
        var species = document.getElementById("species").innerText
        var method = document.getElementById("method").value
        console.log(method)
         
        var logCPM = document.getElementById("logCPM").checked
        var log = document.getElementById("log").checked
        var q = document.getElementById("q").checked
        var z = document.getElementById("z").checked
        var norms = {'logCPM': logCPM, 'log': log, 'q': q, 'z': z}
        console.log(norms)


        var gsedata = JSON.stringify({'gse': gse, 'control': control_condition, 'perturb': perturb_condition, 'method': method, 'species': species, 'norms': norms});
        console.log(gsedata)
        $.ajax({
            url: "/dgeapi",
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
            var api = currURL.filter(x => !x.includes('GSE')).join('/') + '/singlegene/'
            var name = `${control_condition}-vs-${perturb_condition}-${method}`
            var pvals;
            if (method === 'limma') {
                tabletext += "<th></th><th>logFC</th><th>AvgExpr</th><th>t</th><th>P Value</th><th>Adj. P Value</th><th>B</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a href='${api}${vals[0]}' target='_blank'>`  + vals[0] + "</a></td><td>" + vals[1]+"</td><td>" + vals[2] +"</td><td>"+ vals[3] + "</td><td>"+ vals[4] + "</td><td>"+ vals[5] + "</td><td>"+ vals[6] + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'desc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });
                pvals = table.column(4).data()

            } else if (method === 'edgeR') {
                tabletext += "<th></th><th>logFC</th><th>logCPM</th><th>PValue</th><th>FDR</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a href='${api}${vals[0]}' target='_blank'>`  +vals[0] + "</a></td><td>" + vals[1]+"</td><td>" + vals[2] +"</td><td>"+ vals[3] + "</td><td>"+ vals[4] + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext 
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'desc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });
                pvals = table.column(3).data()

            } else if (method === 'DESeq2') {
                tabletext += "<th></th><th>baseMean</th><th>log2FC</th><th>lfcSE</th><th>stat</th><th>P-value</th><th>Adj. P-value</th></tr><tbody>"
                rows.forEach(function(row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a href='${api}${vals[0]}' target='_blank'>` +vals[0] +"</a></td><td>" + vals[1]+"</td><td>" + vals[2] +"</td><td>"+ vals[3] + "</td><td>"+ vals[4] + "</td><td>" + vals[5] + "</td><td>" + vals[6] + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[2, 'desc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', {extend: 'csv', title: name}
                    ]
                });    
                pvals = table.column(5).data()
            }
            var genes = table.column(0).data()

            
            genes = genes.map(x => x.replace(/<\/?[^>]+(>|$)/g, ""))

            var genelist_buttons = 
            `<div class="row justify-content-center mx-auto text-center">
            <div class="mt-3 h7">Submit the </div>
            <select id='dir' class='dirpicker m-2'>
                <option value="top"selected>top</option>
                <option value="bot">bottom</option>
            </select>
            <input id='numgenes' type='number' step='1' value='100' pattern='[0-9]' min='1' class='m-2' style='width: 60px;'/>
            <div class="mt-3 h7">differentially expressed genes with a p-value less than</div>
            <input id='signifigance' type='number' step='.001' value='.05' max='1' class='m-2' style='width: 60px;'/>
            <div class="mt-3 h7">to</div>
            </div>
            <div class="row justify-content-center mx-auto text-center">
            <button class="btn btn-primary btn-group-sm m-2" onclick="submit_geneset('${genes.join(',')}', '${pvals.join(',')}')">Gene Set Queries</button>
            <div class="mt-3 h7">or to</div>
            <button class="btn btn-primary btn-group-sm m-2" onclick="submit_geneset_home('${genes.join(',')}', '${pvals.join(',')}', '${name}')">Diabetes Gene Set Library</button>
            </div>
            `
                
            document.getElementById("geneset-buttons").innerHTML = (clear_button + genelist_buttons)
        });
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


})

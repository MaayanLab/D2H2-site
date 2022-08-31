
// COUNT GENES IN TEXT BOXS ON geneset page

function geneCount(gene_list, num) {
    const genes = gene_list.toUpperCase().split(/\r?\n/g).filter(Boolean);
    $('span#gene-count' + String(num)).text(genes.length);
}

function on_change(el) {

    for (var i =0; i < el.options.length; i++) {
        document.getElementById(el.options[i].value).style.display = 'none';
    }
    document.getElementById(el.options[el.selectedIndex].value).style.display = 'block'; // Show el

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


    
    for (var i = 1; i < 7; i++) {
        var selectize = $(`#search${i}`)[0].selectize;
        selectize.setValue(gene);
    }
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
                table = $('#table-resources').DataTable();
            });

        });
    }

    if (currURL[3] == 'resources') {createResourcesTable();}


    if (currURL[3] == 'downloads') {createDownloadsTable();}

    if (currURL[3] ==  'singlegene' && currURL.length == 5) {
        var gene = currURL[4]
        fillSingleExample(gene);
    } 


    
/*   $('#search1').on('', function () { 
        console.log($('#search1').val())
        fillSingleExample($('#search1').val())
    }) */



    



    


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
    

    $('.search').each(function() {
        $(this).selectize({
            preload: true,
        valueField: 'gene_symbol',
        labelField: 'gene_symbol',
        searchField: 'gene_symbol',
        maxItems: 1,
        render: {
            option: function (item, escape) {
                return '<div class="pt-2 light">'+item.gene_symbol+'</div>';
            }
        },
        load: function (query, callback) {
            // if (!query.length) return callback();
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
        },
        onDropdownClose: function(value) {
            var gene = this.getValue()
            fillSingleExampleSkip(gene, this.id)
        },
    })
    })





      var $gene_select = $('#gene-select').selectize({
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
                }
            });
        },
        persist: false,
    });

    
      

    // CHANGE LINKS FOR APPYTERS/ARCHS4/GTEx DYNAMICALLY   
    
    $('#appyter-home').click( async function() {  
        if ($("#gene-select").val()) {
            var inputvalue = $("#gene-select").val();
            if (document.getElementById('human').className.includes('btn-primary')) {
                var species = 'Human';
                var arg = 'human_gene';
            } else {
                var species = 'Mouse';
                var arg = 'mouse_gene';
            }
            const formData = new FormData()
            formData.append('species_input', species)
            formData.append(arg, inputvalue)
            
            var res = await fetch("https://appyters.maayanlab.cloud/Gene_Expression_T2D_Signatures/", {
                method: "POST",
                headers: {
                    'Accept': 'application/json',
                },
                body: formData,
            })

            const id = await res.json()

            const final_url = "https://appyters.maayanlab.cloud/Gene_Expression_T2D_Signatures/" + id.session_id
            window.open(final_url, '_blank')

            
        } else {
            window.open('https://appyters.maayanlab.cloud/Gene_Expression_T2D_Signatures/', '_blank')
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
        console.log($("#search1").val())
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
                
                tabletext += `<button class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-start text-left"
                        data-toggle="collapse" data-target="#genesoverlap-${data[k][0]}" aria-expanded="false"
                        aria-controls="genesoverlap-${data[k][0]}">
                        <div class="text">Show Overlapping Genes</div>
                </button>
                    <div class="collapse" id="genesoverlap-${data[k][0]}">
                        ${data[k][5].join(", ")}
                    </div></td>
                `
                tabletext += "<td>" + Number(data[k][6]).toPrecision(4) + "</td></tr>";
            }
            tabletext += "</tbody></table>";


            $(document).ready(function(){
                $('#table-enrichr').DataTable();
                    
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
                $('#table-gwas').DataTable();
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
                    $('#tablecor').DataTable();
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
                $('#table-pheno').DataTable();
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

                $(`.table-enrichr`).DataTable();
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

            gene_symbol = 'A1BG';

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

    $('#studies-table').DataTable({responsive: true});


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

        if (!control_condition || !perturb_condition ) {
            alert("Please select both conditions")
            return;
        }

        document.getElementById("dgea-loading").innerHTML = "<div class='loader justify-content-center'></div>";

        var gse = currURL[3]
        var gsedata = JSON.stringify({'gse': gse, 'control': control_condition, 'perturb': perturb_condition});


        
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

            const final_url = "https://appyters.maayanlab.cloud/Bulk_RNA_seq/" + id.session_id + "/#differential-gene-expression"
            window.open(final_url, '_blank')

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


    $('#human').click(function() {
        $('#human').prop('class', 'species-tab btn btn-primary m-1');
        $('#mouse').prop('class', 'species-tab btn btn-inactive m-1');
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
    })

    $('#mouse').click(function() {
        $('#mouse').prop('class', 'species-tab btn btn-primary m-1');
        $('#human').prop('class', 'species-tab btn btn-inactive m-1');
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
    })
    




})

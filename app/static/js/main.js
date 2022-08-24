
// COUNT GENES IN TEXT BOXS ON geneset page

function geneCount(gene_list, num) {
    const genes = gene_list.toUpperCase().split(/\r?\n/g).filter(Boolean);
    $('span#gene-count' + String(num)).text(genes.length);
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


// FILL IN ALL BOXS WITH GENE ENTERED ON singlegene page

function fillPage(event, ui) {

    var gene = ui.item.label;
    const numSearch = 6
    for (let i = 1; i < (numSearch +1); i++) {
        document.getElementById(`search${i}`).value = gene;
    }
    
}

$(document).ready(function() {

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
            console.log()
        });
    }

    if (currURL[3] == 'resources') {createResourcesTable();}


    if (currURL[3] == 'downloads') {$('#table-downloads').DataTable();}



    
    $('a').each(function(){
        
        if ($(this).prop('href') == window.location.href) {
            console.log(($(this).prop('href')));
            $(this).addClass('active'); 
            $(this).parents('li').addClass('active');
        }
    });
    

    
    

    $('.search').autocomplete({
    source: function (request, response) {
        $.ajax({
            url: "/static/searchdata/allgenes.json",
            dataType: 'json',
            data: request,
            success: function( data ) {
                 {
                    var filtered = data.filter(function (str) { 
                        return str.includes(request.term.toUpperCase())});
                    var aSearch = [];
                    // for each element in the main array ...
                    $(filtered).each(function(iIndex, sElement) {
                        // ... if element starts with input value ...
                        if (sElement.substr(0, request.term.length) == request.term.toUpperCase()) {
                            // ... add element
                            aSearch.push(sElement);
                        }
                    });

                    response(aSearch.splice(0,50))
                };
            }
        }); 
       },  
       minLength: 1,
       scroll:true,
       max:5,
       select: fillPage
      });

      $('.search-home').autocomplete({
        source: function (request, response) {
            if (document.getElementById('human').className.includes('btn-primary')) {
                var url =  "/static/searchdata/allgenes.json";
            } else {
                var url =  "/static/searchdata/mouse_genes.json";
            }
            

            $.ajax({
                url: url,
                dataType: 'json',
                data: request,
                success: function( data ) {
                     {
                        var filtered = data.filter(function (str) { 
                            return str.includes(request.term.toUpperCase())});
                        var aSearch = [];
                        // for each element in the main array ...
                        $(filtered).each(function(iIndex, sElement) {
                            // ... if element starts with input value ...
                            if (sElement.substr(0, request.term.length) == request.term.toUpperCase()) {
                                // ... add element
                                aSearch.push(sElement);
                            }
                        });
    
                        response(aSearch.splice(0,50))
                    };
                }
            }); 
           },  
           minLength: 1,
           scroll:true,
           max:5,
    });


      // GET CURRENT URL FOR STUDY SPECIFIC AUTOCOMPLETE

      // DO STUDY SPECIFIC AUTOCOMPLETE

      $('.studysearch').autocomplete({
        source: function (request, response) {
            $.ajax({
                url: "/static/data/" + document.getElementById("species").innerText + "/" + currURL[3] +"/" + "genes.json",
                dataType: 'json',
                data: request,
                success: function( data ) {
                     {
                        var filtered = data.filter(function (str) { 
                            return str.includes(request.term.toUpperCase())});
                        var aSearch = [];
                        // for each element in the main array ...
                        $(filtered).each(function(iIndex, sElement) {
                            // ... if element starts with input value ...
                            if (sElement.substr(0, request.term.length) == request.term.toUpperCase()) {
                                // ... add element
                                aSearch.push(sElement);
                            }
                        });
    
                        response(aSearch.splice(0,50))
                    };
                }
            }); 
           },  
           minLength: 1,
           scroll:true,
           max:5
        });
      
      

    // CHANGE LINKS FOR APPYTERS/ARCHS4/GTEx DYNAMICALLY   
    
    $('#appyter-home').click( async function() {  
        if ($("#search-home").val()) {
            var inputvalue = $("#search-home").val();
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
        if ($("#search1").val()) {
            var inputvalue = $("#search1").val();
            $('#appyter-url1').prop('href', "https://appyters.maayanlab.cloud/Gene_Expression_by_Tissue/#/?args.gene=" +inputvalue + "&submit");
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

    
    // QUERY GWAS AND PRODUCE TABLE

    $('#gwas-query').click(function() {  
        var inputvalue = $("#search4").val();
        if (!inputvalue) {
            $("#gwas-res").html("");
            return;
        }


        $("#gwas-res").addClass('loader');

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

            $("#gwas-res").removeClass('loader');
            document.getElementById("gwas-res").innerHTML = tabletext + clear_button;
        });
    });


    // QUERY ARCHS4 AND CREATE TABLE

    function loadCorrelation(gene){
        $("archs4-res").html("");

        document.getElementById("archs4-res").innerHTML = "<div class='loader justify-content-center mx-auto'> </div>";
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
        $("#komp-res").innerHTML = "<div class='loader justify-content-center'> </div>";

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
        $("#tf-res").innerHTML = "<div class='loader justify-content-center'> </div>";

        $.ajax({
            url: "/gettfs",
            type: "POST",
            data: {gene:inputvalue}
        }).done(function(response) {

            console.log(response)

            const data = response['data'];



            if (data.length === 0) {
                document.getElementById("tf-res").innerHTML = "<p class='text-center'> No data found </p>";
                return;
            }

            var res_html = `<table id='table-enrichr' class='styled-table'><thead><tr><th>Library</th><th></th></tr><tbody>`

            for (var i = 0; i < data.length; i++) {
                var lib = data[i]['name'];
                var sentence = data[i]['format'];

                

                for (j = 0; j < data[i]['tfs'].length; j++) {
                    var tf = data[i]['tfs'][j]
             
                    tf_sentence = sentence.replace('{0}', inputvalue).replace('{1}', tf)
                    
                    res_html += `<tr><td>${lib}</td><td> ${tf_sentence} </td></tr>`

                }
                

            }
            res_html += `</tbody></table>`

            

            $(document).ready(function(){

                $(`#table-enrichr`).DataTable();
            });
            
        

            document.getElementById("tf-res").innerHTML = res_html;

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
    $('#generate-plot').on('click', function(evt) {
        boxplot();
    })
    
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
        console.log(inputvalue)

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

    // PRODUCE INPUT HTML FOR SINGLE GENE SET
    
    function getSingleEntry(num, genecount) {
        var single_entry = `<div class="col-6 text-right mr-3">
            <textarea name="list" rows="8" id="text-area${num}" placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box" onkeyup="geneCount($(this).val(), ${genecount})" onchange="geneCount($(this).val(), ${genecount})" onfocus="geneCount($(this).val(), ${genecount})"></textarea>
            <div class="mt-1">
                <span id="gene-count${genecount}"> 0 </span> gene(s) entered
            </div>
            </div>
            <div class="col-3 justify-content-around">
            <p>
            File formats excepted: csv, tsv, txt file with Entrez gene symbols on each line
            </p>
            <form action="/action_page.php">
        <input type="file" id="gene-file${num}" name="filename">
        </form>
        
        </div>`;
        return single_entry;
    }

    // PRODUCE INPUT HTML FOR MULTIPLE GENE SETS

    function getMultipleEntries(num, genecount) {
        var multiple_entries = `<div class="col-3 text-right">
            <textarea name="list" rows="8" id="text-area${num}-up" placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box" onkeyup="geneCount($(this).val(), ${genecount})" onchange="geneCount($(this).val(),  ${genecount})" onfocus="geneCount($(this).val(),  ${genecount})"></textarea>
            <div class="mt-1">
                <span id="gene-count${genecount}"> 0 </span> UP gene(s) entered
            </div>
            </div>
            <div class="col-3">
            <p>
                File formats excepted: csv, tsv, txt file with Entrez gene symbols on each line
            </p>
            <form action="/action_page.php">
                <input type="file" id="gene-file${num}-up" name="filename">
            </form>
            
            </div>
            <div class="col-3 text-right">
            <textarea name="list" rows="8" id="text-area${num}-down" placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box" onkeyup="geneCount($(this).val(), ${genecount + 1})" onchange="geneCount($(this).val(), ${genecount + 1})" onfocus="geneCount($(this).val(), ${genecount + 1})"></textarea>
            <div class="mt-1">
                <span id="gene-count${genecount + 1}"> 0 </span> DOWN gene(s) entered
            </div>
            </div>
            <div class="col-3">
            <p>
                File formats excepted: csv, tsv, txt file with Entrez gene symbols on each line
            </p>
            <form action="/action_page.php">
                <input type="file" id="gene-file${num}-down" name="filename">
            </form>
            
            </div>`;
        return multiple_entries;
    }

    // FOR KEA3/CHEA3 INPUT SWITCHING  (between multiple and single gene set)

    /* $('#ea3_entries_button').click( function() {  
        
        var button = document.getElementById("ea3_entries_button");

        if (button.textContent === "Use Up/Down Gene Sets"){
            document.getElementById("ea3_entries").innerHTML = getMultipleEntries(2, 2);
            button.textContent = "Use Single Gene Set"
        } else {
            document.getElementById("ea3_entries").innerHTML = getSingleEntry(2, 2);
            button.textContent = "Use Up/Down Gene Sets"
        }
        
    }); */

    // SWITCH INPUT FOR SIGCOM LINCS

    $('#sigcom_entries_button').click( function() {  
        
        var button = document.getElementById("sigcom_entries_button");

        if (button.textContent === "Use Up/Down Gene Sets"){
            document.getElementById("sigcom_entries").innerHTML = getMultipleEntries(3, 4);
            button.textContent = "Use Single Gene Set"
        } else {
            document.getElementById("sigcom_entries").innerHTML = getSingleEntry(3, 4);
            button.textContent = "Use Up/Down Gene Sets"
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

        var gse = currURL[3]
        var gsedata = JSON.stringify({'gse': gse, 'control': control_condition, 'perturb': perturb_condition});
        
        $.ajax({
            url: "/api/data",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done( async function(response) {
            console.log(response)
            console.log(response['meta'])

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

            console.log(res)
            const id = await res.json()

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
        $('#human').prop('class', 'species-tab btn btn-primary m-1')
        $('#mouse').prop('class', 'species-tab btn btn-inactive m-1')
    })

    $('#mouse').click(function() {
        $('#mouse').prop('class', 'species-tab btn btn-primary m-1')
        $('#human').prop('class', 'species-tab btn btn-inactive m-1')
    })
    




})

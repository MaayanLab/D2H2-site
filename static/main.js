


$(document).ready(function() {
            

    $('.search').autocomplete({
    source: function (request, response) {
        $.ajax({
            url: "/static/genes.json",
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

                    response(aSearch)
                };
            }
        }); 
       },  
       minLength: 1,
       scroll:true,
       max:5
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


    $('#gwas-query').click(function() {  
        var inputvalue = $("#search4").val();
        if (!inputvalue) {
            $("#gwas-res").html("");
            return;
        }

        $.ajax({
            url: "/getgwas",
            type: "POST",
            data: {gene:inputvalue}
        }).done(function(response) {

            const data = response['GWAS_Catalog']
 
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


    function loadCorrelation(gene){
        $("archs4-res").html("");

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

    // Is there a knockout mouse ==> query komp

    function hideEvidence() { 
        var table = $('#table-pheno')
        console.log(table)
 
        // Get the column API object
        var column = table.column$(3);
 
        // Toggle the visibility
        column.visible(!column.visible());
    };

    $('#komp-query').click(function() {  
        var inputvalue = $("#search6").val();
        if (!inputvalue) {
            $("#komp-res").html("");
            return;
        }

        $.ajax({
            url: "/getkomp",
            type: "POST",
            data: {gene:inputvalue}
        }).done(function(response) {

                //row["subject.primaryIdentifier"], row["subject.symbol"], \
                //row["subject.sequenceOntologyTerm.name"], row["ontologyTerm.identifier"], \
                //row["ontologyTerm.name"], row["evidence.publications.pubMedId"], \
                //row["evidence.comments.type"], row["evidence.comments.description"])

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


    


});


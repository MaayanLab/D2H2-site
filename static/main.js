


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

                    response(aSearch.splice(0, 6))
                };
            }
        }); 
       },  
       minLength: 1
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

            const data = response['GWAS_Catalog_2019']
 
            if (data.length === 0) {
                $("#gwas-res").html("<p class='text-center'> No data found </p>");
                return;
            }

            
            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#gwas-res\").html(\"\");'> Clear Results </button> </a>"
            const header = ["Name", "p-value", "Adjusted p-value", "Odds ratio", "Combined score"]
            const cols = [1, 2, 6, 3, 4]

            var html = "<table class='styled-table' id='table-gwas'><thead><tr>"

            for (val in header) {
                html += ('<th>' + header[val] + '</th>');
            }

            html += ("</tr><tbody>");

            for (i in data) {
                html += ('<tr>');
                for (j in cols) {

                    var val = data[i][cols[j]];
                    if (cols[j] === 2 || cols[j] === 6) {
                        val = Number(val).toPrecision(4)
                    }

                    if (cols[j] === 3 || cols[j] === 4) {
                        val = Number(val).toFixed(2)
                    }

                    html += '<td>' + val + '</td>';
                }
                html += ('</tr>');
            }
            html += ("</tbody></table>");

            $(document).ready(function(){
                $('#table-gwas').DataTable();
            });

            document.getElementById("gwas-res").innerHTML = html + clear_button;
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
                    tabletext += "<tr><td>"+k+"</td><td><a href=\"https://maayanlab.cloud/archs4/gene/"+genesym[k]+"\" target=\"_blank\">"+genesym[k]+"</a></td><td>"+correlation[k]+"</td></tr>";
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

            const num_pheno = response['response']['numFound']
            const data = response['response']['docs']

            if (num_pheno === 0) {
                $("#komp-res").html("<p class='text-center'> No data found </p>");
                return;
            }
            const clear_button = "<a> <button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='$(\"#komp-res\").html(\"\");'> Clear Results </button> </a>"
            var tabletext = "<table id='table-pheno' class='styled-table'><thead><tr><th>Phenotype</th><th>System</th><th>p-value</th><th>Source</th></tr><tbody>";
            for (var k = 0; k < num_pheno; k++) {
                tabletext += "<tr><td><a href='https://www.mousephenotype.org/data/genes/" +data[k]['marker_accession_id'] + "' target='_blank'>"+data[k]['mp_term_name']+"</td><td>"+data[k]['top_level_mp_term_name'].join(" ")+"</td><td>"+data[k]['p_value']+"</td><td>"+data[k]['pipeline_name']+"</td></a></tr>";
            }
            tabletext += "</tbody></table>";

            $(document).ready(function(){
                $('#table-pheno').DataTable();
            });

            document.getElementById("komp-res").innerHTML = tabletext + clear_button;

        });
    });


    


});


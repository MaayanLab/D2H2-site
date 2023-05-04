import { human_list, mouse_list, processes } from "./constants.js";
import { initialize_search, gen_table, gene_signatures, generanger_plot, geo_reverse, single_gene_perturbations, l1000_reverse, query_gwas, query_enrichr_tfs, loadCorrelation, query_komp } from './single-gene-queries.js';
import { geneset_signatures, geneset_enrichment, geneset_kea3, geneset_chea3, geneset_sigcomlincs } from './geneset-queries.js';

export function runFindQuery(q) {
    clear_home()
    var data = JSON.stringify({ 'query': q })
    $.ajax({
        url: "api/query_gpt",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: data
    }).done(async function (response) {
        if (response['response'] == 1) {
            prompt('Error analyzing query... Please try again')
            return;
        }
        var gene = '';
        var species = '';
        if (response['input'] == '[Gene]') {
            var check_list = await human_list
            var check_list_mouse = await mouse_list
            q.split(' ').forEach((w) => {
                var w_clean = w.replace("?", '').replace('.', '')
                if (check_list.includes(w_clean)) {
                    gene = w_clean;
                    species = 'Human'
                }
                else if (check_list_mouse.includes(w_clean)) {
                    gene = w_clean;
                    species = 'Mouse'
                }
            })

            if (gene == '') {
                initialize_search('');
                document.getElementById("output").innerText = response['output'];
                document.getElementById("selector").style.display = "flex";
                document.getElementById("submit_single_gene").style.display = "flex";
            } else {
                document.getElementById('loading').innerHTML = "<div class='loadingspinner'><div id='square1'></div><div id='square2'></div><div id='square3'></div><div id='square4'></div><div id='square5'></div></div>"
                run_process_gene(gene, species, response['output']).then(() => document.getElementById('loading').innerHTML = "")
            }
        } else if (response['input'] == '[GeneSet]') {
            //window.open(window.location.href + 'geneset#' + response['output'], "_self")
            document.getElementById("output").innerText = response['output'];
            document.getElementById("enter-geneset").style.display = "flex";
            if (response['output'] == '[L1000]') document.getElementById("toggle").style.display = 'flex';
            document.getElementById("submit_gene_set").style.display = 'flex';

        } else document.getElementById('loading').innerHTML = "";
    });
}



export async function submit_single_gene() {
    var gene = $("#search").val()
    if (!gene) {
        alert('Enter gene symbol')
        return;
    } else {
        var check_list = await human_list
        var check_list_mouse = await mouse_list
        var species;
        if (check_list.includes(gene)) {
            species = 'Human'
        }
        else if (check_list_mouse.includes(gene)) {
            species = 'Mouse'
        }
        var output = document.getElementById("output").innerText;
        clear_home();
        document.getElementById('loading').innerHTML = "<div class='loadingspinner'><div id='square1'></div><div id='square2'></div><div id='square3'></div><div id='square4'></div><div id='square5'></div></div>";
        run_process_gene(gene, species, output).then(() => document.getElementById('loading').innerHTML = "");
    }
}

export async function run_process_gene(gene, species, output) {
    var process_eval = await processes;
    var process = process_eval["[Gene]"][output];
    var args = process.args.join(",");
    if (process.questions.length > 0) {

    }
    document.getElementById("description").innerText = process.description;
    args = args.replace("geneSymbol", gene);
    args = args.replace("species", species);
    args = args.split(',')
    var args_string = "";
    for (let i = 0; i < args.length; i++) {
        if (i == 0) args_string += `'${args[i]}'`;
        else args_string += `,'${args[i]}'`;
    }
    var process_eval = `${process.process}(${args_string})`;
    await eval(process_eval);
    return
}

export async function submit_gene_set() {
    var geneset = '';
    var geneset_up = '';
    var geneset_down = '';

    if (document.getElementById("enter-geneset").style.display == "flex") {
        var inputvalue = document.getElementById("text-area").value;
        var file = document.getElementById("gene-file").value;
        var section = "gene-file";
        geneset;

        if (inputvalue) {
            geneset = inputvalue.split("\n").join(',')
        }

        if (file) {
            geneset = await loadFileAsText(section, ",");
        }
    } else {
        var inputvalue_up = document.getElementById("text-area-up").value;
        var inputvalue_down = document.getElementById("text-area-down").value;
        var file_up = document.getElementById("gene-file-up").value;
        var file_down = document.getElementById("gene-file-down").value;
        var section_up = "gene-file-up";
        var section_down = "gene-file-down";

        if (inputvalue_up) {
            geneset_up = inputvalue_up.split("\n").join(',')
        }
        if (file_up) {
            geneset_up = await loadFileAsText(section_up, ",");
        }

        if (inputvalue_down) {
            geneset_down = inputvalue_down.split("\n").join(',')
        }

        if (file_down) {
            geneset_down = await loadFileAsText(section_down, ",");
        }
    }
    var output = document.getElementById("output").innerText;
    clear_home();
    document.getElementById('loading').innerHTML = "<div class='loadingspinner'><div id='square1'></div><div id='square2'></div><div id='square3'></div><div id='square4'></div><div id='square5'></div></div>";
    run_process_geneset(geneset, geneset_up, geneset_down, output).then(() => document.getElementById('loading').innerHTML = "");
}

export async function run_process_geneset(geneset, geneset_up, geneset_down, output) {
    var process_eval = await processes;
    var process = process_eval["[GeneSet]"][output];
    var args = process.args.join("\t");
    if (process.questions.length > 0) {

    }
    document.getElementById("description").innerText = process.description;
    args = args.replace("geneset", geneset);
    args = args.replace("up", geneset_up);
    args = args.replace("down", geneset_down);
    args = args.split('\t')
    var args_string = "";
    for (let i = 0; i < args.length; i++) {
        if (i == 0) args_string += `'${args[i]}'`;
        else args_string += `,'${args[i]}'`;
    }
    var process_eval = `${process.process}(${args_string})`;
    console.log(process_eval)
    await eval(process_eval);
    return
}
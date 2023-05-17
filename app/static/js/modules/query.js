import { human_list, mouse_list, processes, chatN, chatNresult, loading, geneset_entries} from "./constants.js";
import { gen_table, gene_signatures, generanger_plot, geo_reverse, single_gene_perturbations, l1000_reverse, query_gwas, query_enrichr_tfs, loadCorrelation, query_komp } from './single-gene-queries.js';
import { geneset_signatures, geneset_enrichment, geneset_kea3, geneset_chea3, geneset_sigcomlincs } from './geneset-queries.js';


export async function select_option(q, options) {
    var data = JSON.stringify({ 'response': q , 'options': options})
    var response = await $.ajax({
        url: "api/query_options",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: data
    })
    return response['option']
}

export async function parse_gene(q) {
    var gene = '';
    var species = '';
    var check_list = await human_list
    var check_list_mouse = await mouse_list
    q.split(' ').forEach((w) => {
        var w_clean = w.replace("?", '').replace('.', '')
        if (check_list.includes(w_clean)) {
            gene = w_clean;
            species = 'Human';
        }
        else if (check_list_mouse.includes(w_clean)) {
            gene = w_clean;
            species = 'Mouse';
        }
    })
    return [gene, species]
}


export async function runFindQuery(q) {
    var data = JSON.stringify({ 'query': q })
    var response = await $.ajax({
        url: "api/query_gpt",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: data
    })
    if (response['response'] == 1) {
        return response;
    }
    if (response['input'] == '[Gene]') {
        const res = await parse_gene(q);
        response['geneSymbol'] = res[0]
        response['species'] = res[1]
        return response

    } else if (response['input'] == '[GeneSet]') {
        response['geneset'] = '';
        return response

    } else document.getElementById('loading').innerHTML = "";
}


export async function run_process_gene(process_info_copy, chat_num) {
    var chat_num = chat_num;
    var process_eval = await processes;
    var process = process_eval[process_info_copy.input][process_info_copy.output];
    var args = process.args.map((x) => x);
    if (process.questions.length > 0) {

    }
    document.getElementById("chat-bubbles-section").appendChild(chatN('start', chat_num, '#d3d3d3', process.text))
    $(`#chat-${chat_num}`).fadeIn(2000, async () => {
        chat_num++;
        const placeholder = document.createElement("div");
        placeholder.innerHTML = `<div id='loading${chat_num}'>${loading}</div>`;
        const loadingNode = placeholder.firstElementChild;
        document.getElementById("chat-bubbles-section").appendChild(loadingNode)
        document.getElementById("chat-bubbles-section").appendChild(chatNresult('start', chat_num, '#d3d3d3', "result" + chat_num))
        
        process_info_copy['result'] = "result" + chat_num
        for (let j = 0; j < args.length; j ++) {
            args[j] = process_info_copy[args[j]];
        }
        var args_string = "";
        for (let i = 0; i < args.length; i++) {
            if (i == 0) args_string += `'${args[i]}'`;
            else args_string += `,'${args[i]}'`;
        }
        var process_eval = `${process.process}(${args_string})`;
        console.log(process_eval)
        await eval(process_eval);
        document.getElementById('loading' + chat_num).style.display = 'none';
        await $(`#chat-${chat_num}`).fadeIn(2000);
    
        return;
    })
    
}

export const human_list = fetch(
    "static/data/allgenes-comb.json"
).then(data => data.json());

export const mouse_list = fetch(
    "static/data/mouse_genes.json"
).then(data => data.json());

export const processes = fetch(
    "static/data/processes.json"
).then(data => data.json());

export const loading = "<div class='loadingspinner'><div id='square1'></div><div id='square2'></div><div id='square3'></div><div id='square4'></div><div id='square5'></div></div>";




export function chatNresult(side, chatNum, color, resultid) {
    var icon;
    var padding;
    if (side == 'start') {
        icon = "static/img/d2h2_chat_bot.png";
        padding = "pl-2";
    }
    else {
        icon = "static/img/user_icon.png";
        padding = "pr-2";
    }
    const chathtml = `
    <div class="chat chat-${side} mb-1" id="chat-${chatNum}" style="display: none; opacity: 1 !important;"> 
        <div class="chat-image avatar ${padding}" style="width: 50px">
            <img src="${icon}"/>
        </div>
        <div class="chat-header">
            D2H2
        </div>
        <div class="chat-bubble chat-bubble-info" id="${resultid}" style="background-color: ${color}; overflow: auto;">
        </div>
    </div>
    `   
    const placeholder = document.createElement("div");
    placeholder.innerHTML = chathtml;
    const node = placeholder.firstElementChild;
    return node
}


export function geneset_entries(chatnum) { 
    return `
    <div class="row justify-content-center" id="enter-geneset${chatnum}" style="display: flex;">
        <div class="col-sm-12 col-md-6 col-lg-3 text-center justify-content-end mr-4">
            <textarea class="input-form" name="list" rows="8" id="text-area${chatnum}"
                placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                onkeyup="geneCount($(this).val(), ${chatnum})" onchange="geneCount($(this).val(), ${chatnum})"
                onfocus="geneCount($(this).val(), ${chatnum})"></textarea>
            <div class="mt-1">
                <span id="gene-count${chatnum}"> 0 </span> gene(s) entered
            </div>
            <div class="text-text-center">
                <a onclick="fillSet('text-area${chatnum}', 'desc', ${chatnum}); "
                    style="color: rgb(10, 13, 149)">Try an example set</a>
            </div>
        </div>
        <div class="col-sm-12 col-md-6 col-lg-3 justify-content-start text-center ml-4"
            style="overflow-y: visible !important;">
            <p>File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line</p>
            <form action="/action_page.php">
                <input type="file" id="gene-file${chatnum}" name="filename">
            </form>
            <label class="mt-1 mr-1" for="desc">Enter gene set description (optional): </label>
            <input class="form-control-md" id="desc" type="text">
        </div>
    </div>`
}

export function geneset_up_down_entries(chatnum) { 
    return `
    <div class="row justify-content-center" id="enter-geneset${chatnum}" style="display: none;">
        <div class="col-sm-12 col-md-6 col-lg-3 text-center justify-content-end mr-4">
            <textarea class="input-form" name="list" rows="8" id="text-area${chatnum}"
                placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                onkeyup="geneCount($(this).val(), ${chatnum})" onchange="geneCount($(this).val(), ${chatnum})"
                onfocus="geneCount($(this).val(), ${chatnum})"></textarea>
            <div class="mt-1">
                <span id="gene-count${chatnum}"> 0 </span> gene(s) entered
            </div>
            <div class="text-text-center">
                <a onclick="fillSet('text-area${chatnum}', 'desc', ${chatnum}); "
                    style="color: rgb(10, 13, 149)">Try an example gene set</a>
            </div>
        </div>
        <div class="col-sm-12 col-md-6 col-lg-3 justify-content-start text-center ml-4"
            style="overflow-y: visible !important;">
            <p>File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line</p>
            <form action="/action_page.php">
                <input type="file" id="gene-file${chatnum}" name="filename">
            </form>
            <label class="mt-1 mr-1" for="desc">Enter gene set description (optional): </label>
            <input class="form-control-md" id="desc" type="text">
        </div>
    </div>
    <div class="flex-row justify-content-center" id="enter-geneset-up-down${chatnum}"
        style="flex-wrap: wrap; display: flex;">
        <div class="col">
            <div class="col text-center mr-3">
                <textarea class="input-form" name="list" rows="8" id="text-area-up${chatnum}"
                    placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                    onkeyup="geneCount($(this).val(), ${chatnum + 1})" onchange="geneCount($(this).val(), ${chatnum + 1})"
                    onfocus="geneCount($(this).val(), ${chatnum + 1})"></textarea>
                <div class="mt-1">
                    <span id="gene-count${chatnum + 1}"> 0 </span> UP gene(s) entered
                </div>
                <div class="text-center">
                    <a onclick="fillSet('text-area-up${chatnum}', '', ${chatnum + 1})" style="color: rgb(10, 13, 149)">
                        Try an example gene set</a>
                </div>
            </div>
            <div class="col justify-content-center text-center m-2" style="overflow-y: visible !important;">
                <p>
                    File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line
                </p>
                <form action="/action_page.php">
                    <input type="file" id="gene-file-up${chatnum}" name="filename">
                </form>
            </div>
        </div>
        <div class="col">
            <div class="col text-center mr-3">
                <textarea class="input-form" name="list" rows="8" id="text-area-down${chatnum}"
                    placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                    onkeyup="geneCount($(this).val(), ${chatnum + 2})" onchange="geneCount($(this).val(), ${chatnum + 2})"
                    onfocus="geneCount($(this).val(), ${chatnum + 2})"></textarea>
                <div class="mt-1">
                    <span id="gene-count${chatnum + 2}"> 0 </span> DOWN gene(s) entered
                </div>
                <div class="text-center">
                    <a onclick="fillSet2('text-area-down${chatnum}', '', ${chatnum + 2})" style="color: rgb(10, 13, 149)">Try an
                        example
                        gene set</a>
                </div>
            </div>
            <div class="col justify-content-center text-center m-2" style="overflow-y: visible !important;">
                <p>
                    File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line
                </p>
                <form action="/action_page.php">
                    <input type="file" id="gene-file-down${chatnum}" name="filename">
                </form>
            </div>
        </div>
    </div>`
}




export function geneset_buttons(chatnum){
    return `
<div class="text-center justify-content-center row mb-5">
    <button class='btn btn-primary btn-group-sm text-center' id="submit_gene_set${chatnum}">Submit</button>
</div>`
}

export function geneset_buttons_up_down(chatnum){ 
    return `
<div class="text-center justify-content-center row mb-5">
    <button class='btn btn-primary btn-group-sm text-center mr-3' id="toggle${chatnum}" onclick="up_down_toggle('${chatnum}')">Singe Gene Set</button>
    <button class='btn btn-primary btn-group-sm text-center' id="submit_gene_set${chatnum}">Submit</button>
</div>`
}



export const human_list = fetch(
    "static/searchdata/allgenes-comb.json"
).then(data => data.json());

export const mouse_list = fetch(
    "static/searchdata/mouse_genes.json"
).then(data => data.json());

export const processes = fetch(
    "static/searchdata/processes.json"
).then(data => data.json());

export const loading = "<div class='loadingspinner'><div id='square1'></div><div id='square2'></div><div id='square3'></div><div id='square4'></div><div id='square5'></div></div>";

export function chatN(side, chatNum, color, content) {
    var icon;
    var padding;
    var name;
    if (side == 'start') {
        icon = "static/img/d2h2_chat_bot.png";
        padding = "pl-2";
        name = 'D2H2'
    }
    else {
        icon = "static/img/user_icon.png";
        padding = "pr-2";
        name = ''
    }
    const chathtml = `
    <div class="chat chat-${side} mb-1" id="chat-${chatNum}" style="display: none;"> 
        <div class="chat-image avatar ${padding}" style="width: 50px">
            <img src="${icon}"/>
        </div>
        <div class="chat-header">
            ${name}
        </div>
        <div class="chat-bubble chat-bubble-info" style="background-color: ${color};">
            ${content}
        </div>
    </div>
    `   
    const placeholder = document.createElement("div");
    placeholder.innerHTML = chathtml;
    const node = placeholder.firstElementChild;
    return node
}


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
    <div class="chat chat-${side} mb-1" id="chat-${chatNum}" style="display: none;"> 
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


export const geneset_entries = `<div class="row justify-content-center" id="enter-geneset" style="display: flex;">
                                    <div class="col-sm-12 col-md-6 col-lg-3 text-center justify-content-end mr-4">
                                        <textarea class="input-form" name="list" rows="8" id="text-area"
                                            placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                                            onkeyup="geneCount($(this).val(), 0)" onchange="geneCount($(this).val(), 0)"
                                            onfocus="geneCount($(this).val(), 0)"></textarea>
                                        <div class="mt-1">
                                            <span id="gene-count0"> 0 </span> gene(s) entered
                                        </div>
                                        <div class="text-text-center">
                                            <a id="fill-text-area"
                                                style="color: rgb(10, 13, 149)">Try an example gene set</a>
                                        </div>
                                    </div>
                                    <div class="col-sm-12 col-md-6 col-lg-3 justify-content-start text-center ml-4"
                                        style="overflow-y: visible !important;">
                                        <p>File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line</p>
                                        <form action="/action_page.php">
                                            <input type="file" id="gene-file" name="filename">
                                        </form>
                                        <label class="mt-1 mr-1" for="desc">Enter gene set description (optional): </label>
                                        <input class="form-control-md" id="desc" type="text">
                                    </div>
                                    </div>
                                    <div class="flex-row justify-content-center" id="enter-geneset-up-down"
                                    style="flex-wrap: wrap; display: none;">
                                    <div class="col">
                                        <div class="col text-center mr-3">
                                            <textarea class="input-form" name="list" rows="8" id="text-area-up"
                                                placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                                                onkeyup="geneCount($(this).val(), 1)" onchange="geneCount($(this).val(), 1)"
                                                onfocus="geneCount($(this).val(), 1)"></textarea>
                                            <div class="mt-1">
                                                <span id="gene-count1"> 0 </span> UP gene(s) entered
                                            </div>
                                            <div class="text-center">
                                                <a onclick="fillSet('text-area-up', '', 'up')" style="color: rgb(10, 13, 149)">
                                                    Try an example gene set</a>
                                            </div>
                                        </div>
                                        <div class="col justify-content-center text-center m-2" style="overflow-y: visible !important;">
                                            <p>
                                                File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line
                                            </p>
                                            <form action="/action_page.php">
                                                <input type="file" id="gene-file-up" name="filename">
                                            </form>
                                        </div>
                                    </div>
                                    <div class="col">
                                        <div class="col text-center mr-3">
                                            <textarea class="input-form" name="list" rows="8" id="text-area-down"
                                                placeholder="Paste a set of valid Entrez gene symbols (e.g. STAT3) on each row in the text-box"
                                                onkeyup="geneCount($(this).val(), 2)" onchange="geneCount($(this).val(), 2)"
                                                onfocus="geneCount($(this).val(), 2)"></textarea>
                                            <div class="mt-1">
                                                <span id="gene-count2"> 0 </span> DOWN gene(s) entered
                                            </div>
                                            <div class="text-center">
                                                <a onclick="fillSet2('text-area-down', '', 'up')" style="color: rgb(10, 13, 149)">Try an
                                                    example
                                                    gene set</a>
                                            </div>
                                        </div>
                                        <div class="col justify-content-center text-center m-2" style="overflow-y: visible !important;">
                                            <p>
                                                File formats accepted: csv, tsv, txt file with Entrez gene symbols on each line
                                            </p>
                                            <form action="/action_page.php">
                                                <input type="file" id="gene-file-down" name="filename">
                                            </form>
                                        </div>
                                    </div>
                                </div>`
export const geneset_buttons = `
<div class="text-center justify-content-center row mb-5">
    <button class='btn btn-primary btn-group-sm text-center' id="submit_gene_set">Submit</button>
</div>`

export const geneset_buttons_up_down = `
<div class="text-center justify-content-center row mb-5">
    <button class='btn btn-primary btn-group-sm text-center mr-3' id="toggle" onclick="up_down_toggle()">Up & Down Gene Sets</button>
    <button class='btn btn-primary btn-group-sm text-center' id="submit_gene_set">Submit</button>
</div>`
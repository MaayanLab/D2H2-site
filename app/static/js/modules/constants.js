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
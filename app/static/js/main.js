function up_down_toggle(chatnum) {
    if (document.getElementById(`enter-geneset${chatnum}`).style.display == "flex") {
        document.getElementById(`enter-geneset${chatnum}`).style.display = "none";
        document.getElementById(`enter-geneset-up-down${chatnum}`).style.display = "flex";
        document.getElementById(`toggle${chatnum}`).innerText = "Use Single Gene Set";
    } else {
        document.getElementById(`enter-geneset-up-down${chatnum}`).style.display = "none";
        document.getElementById(`enter-geneset${chatnum}`).style.display = "flex";
        document.getElementById(`toggle${chatnum}`).innerText = "Up & Down Gene Sets";
    }
}

function chatN(side, chatNum, color, content) {
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
        name = 'User'
    }
    const chathtml = `
    <div class="chat chat-${side} mb-1" id="chat-${chatNum}" style="display: none; opacity: 1.0 !important;"> 
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
function make_option(chat_num, option) {
    document.getElementById('chat-bubbles-section').appendChild(chatN('end', chat_num, '#d3d3d3', option))
}



function check_genes_present(genes) {
    if (genes.length == 0) {
        alert('no significantly differentially expressed genes were identified with these thresholds');
        return true;
    } else return false;
}

function filter_genes_sigs(genelist, adjpvals, pvals, logfc) {
    var col = 'pval'
    try {
        col = document.getElementById('col-to-use').value
    } catch {
        col = 'pval'
    }

    genelist = genelist.split(',')
    if (col == 'pval') {
        var sigs = pvals;

    } else var sigs = adjpvals;
    sigs = sigs.split(',').map(function (item) {
        return parseFloat(item);
    });
    logfc = logfc.split(',').map(function (item) {
        return parseFloat(item);
    });

    var signifigance = document.getElementById('signifigance').value

    var dir = document.getElementById('dir').value
    var genes_valid = []
    if (dir == 'up') {
        for (i = 0; i < genelist.length; i++) {
            if (sigs[i] <= signifigance && logfc[i] > 0) {
                genes_valid.push(genelist[i])
            }
        }
    } else {
        for (i = 0; i < genelist.length; i++) {
            if (sigs[i] <= signifigance && logfc[i] < 0) {
                genes_valid.push(genelist[i])
            }
        }
    }
    
    return genes_valid
}


function generate_single_plots() {
    // This function will generate the umap, tsne, and pca plots for each indivdual study for a specific condition
    document.getElementById("umap-plot").innerHTML = "";
    document.getElementById("tsne-plot").innerHTML = "";
    document.getElementById("pca-plot").innerHTML = "";

    // document.getElementById("singleplots-loading").innerHTML = "<div class='loader justify-content-center'></div>";
    $('#singleplots-loading').addClass('loader justify-content-center');
    var gse = document.getElementById("singlegse").innerText
    var species = document.getElementById("species").innerText
    var condition_group = document.getElementById("methodsingle").value

    var gsedata = JSON.stringify({ 'gse': gse, 'species': species, 'conditiongroup': condition_group });
    $.ajax({
        url: "/singleplots",
        contentType: 'application/json',
        type: "POST",
        dataType: 'json',
        data: gsedata,
    }).done(function (response) {
        document.getElementById("singleplots-loading").innerHTML = "";

        $('#singleplots-loading').removeClass('loader justify-content-center');
        var umap_plot = response['umapplot']
        var pca_plot = response['pcaplot']
        var tsne_plot = response['tsneplot']

        window.Bokeh.embed.embed_item(umap_plot)
        window.Bokeh.embed.embed_item(pca_plot)
        window.Bokeh.embed.embed_item(tsne_plot)
    });
    $("#boxplot").attr("data-url-plot", `/api/plot_single/${gse}/${condition_group}`)
}


///////// anitmated number counters /////////////
// How long you want the animation to take, in ms
const animationDuration = 600;
// Calculate how long each ‘frame’ should last if we want to update the animation 60 times per second
const frameDuration = 1000 / 60;
// Use that to calculate how many frames we need to complete the animation
const totalFrames = Math.round(animationDuration / frameDuration);
// An ease-out function that slows the count as it progresses
const easeOutQuad = t => t * (2 - t);

// The animation function, which takes an Element
const animateCountUp = el => {
    let frame = 0;
    const countTo = parseInt(el.innerHTML, 10);
    // Start the animation running 60 times per second
    const counter = setInterval(() => {
        frame++;
        // Calculate our progress as a value between 0 and 1
        // Pass that value to our easing function to get our
        // progress on a curve
        const progress = easeOutQuad(frame / totalFrames);
        // Use the progress value to calculate the current count
        const currentCount = Math.round(countTo * progress);

        // If the current count has changed, update the element
        if (parseInt(el.innerHTML, 10) !== currentCount) {
            el.innerHTML = currentCount;
        }

        // If we’ve reached our last frame, stop the animation
        if (frame === totalFrames) {
            clearInterval(counter);
        }
    }, frameDuration);
};


// COUNT GENES IN TEXT BOXS ON geneset page
function geneCount(gene_list, num) {
    const genes = gene_list.toUpperCase().split(/\r?\n/g).filter(Boolean);
    $('span#gene-count' + String(num)).text(genes.length);
}


function clear_home() {
    try {
    document.getElementById("result").innerHTML = "";
    document.getElementById("selector").style.display = "none";
    document.getElementById("output").innerText = "";
    document.getElementById("description").innerText = "";
    document.getElementById("actions").innerHTML = "";
    document.getElementById("enter-geneset").style.display = "none";
    document.getElementById("enter-geneset-up-down").style.display = "none";
    document.getElementById("toggle").style.display = "none";
    document.getElementById("submit_single_gene").style.display = "none";
    document.getElementById("submit_gene_set").style.display = "none";
    } catch {
    document.getElementById('[Signatures]result').innerHTML = "";
    document.getElementById('[TFs]result').innerHTML = "";
    document.getElementById('[Traits]result').innerHTML = "";
    document.getElementById('[Knockout]result').innerHTML = "";
    document.getElementById('[Correlation]result').innerHTML = "";
    }
}

function submit_geneset(genelist, adjpvals, pvals, logfc) {

    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals, logfc)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value


    if (check_genes_present(genes_valid)) return;


    var genes = genes_valid.splice(0, numgenes).join('&')

    const currURL = window.location.href.split('/')
    localStorage.setItem('geneset', genes)
    var url = currURL.filter(x => !x.includes('GSE')).join('/') + '/geneset'
    window.open(url, '_blank')
}

function setGenes(genes) {
    localStorage.setItem('geneset', genes)
}

function clear_dge() {
    document.getElementById("dge-table-area").innerHTML = ""
    document.getElementById("dge-plot").innerHTML = ""
    document.getElementById("dge-loading").innerHTML = ""
    document.getElementById("geneset-buttons").innerHTML = ""

}
function clear_dge_single() {
    document.getElementById("dge-table-area").innerHTML = ""
    document.getElementById("dge-plot").innerHTML = ""
    document.getElementById("dge-loading").innerHTML = ""
    document.getElementById("geneset-buttons").innerHTML = ""
    document.getElementById("enrichment-area").innerHTML = ""

}

function on_change(el) {

    for (var i = 0; i < el.options.length; i++) {
        document.getElementById(el.options[i].value).style.display = 'none';
    }
    document.getElementById(el.options[el.selectedIndex].value).style.display = 'block'; // Show el

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
            error: (error) => {
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

    var description = options.description || "",
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

function filter_and_submit_to_enrichr(genelist, adjpvals, pvals, logfc, description) {
    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals, logfc)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value
    var genes = genes_valid.splice(0, numgenes).join('\n')

    options = {}
    options.list = genes
    options.description = description
    options.popup = true
    enrich(options)
}

async function filter_and_submit_to_kg(genelist, adjpvals, pvals, logfc, description) {

    var genes_valid = filter_genes_sigs(genelist, adjpvals, pvals, logfc)

    if (check_genes_present(genes_valid)) return;

    var numgenes = document.getElementById('numgenes').value

    var genes_str = genes_valid.splice(0, numgenes).join('\n')

    const formData = new FormData()
    formData.append('list', genes_str)
    formData.append('description', description)

    var res = await fetch("https://maayanlab.cloud/Enrichr/addList", {
        method: "POST",
        headers: {
            'Accept': 'application/json',
        },
        body: formData,
    })

    const result = await res.json()
    const userListId = result['userListId']
    const url = `https://maayanlab.cloud/enrichr-kg?userListId=${userListId}`
    window.open(url, '_blank')
}


function loadFileAsText(section, delim) {

    return new Promise((resolve, reject) => {
        var fileToLoad = document.getElementById(section).files[0];

        var fileReader = new FileReader();
        fileReader.onload = async function (fileLoadedEvent) {
            var textFromFileLoaded = fileLoadedEvent.target.result;
            var genes = textFromFileLoaded.split(/[\n\t,]/).join(delim)
            resolve(genes)
        };

        fileReader.readAsBinaryString(fileToLoad, "UTF-8");
    });
}

function fillSingleExample(gene) {
    $(document).ready(function () {
        document.getElementById("singlegenenav").classList.add('active')
        for (var i = 0; i < 7; i++) {

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
    $('.input-form').each(function () {
        document.getElementById(this.id).value = geneset;
        geneCount(geneset, this.id[this.id.length - 1])
    })
}


function fillSet(id, descid, count_id) {
    $.ajax({
        url: "getexample",
        type: "POST",
        data: {},
        dataType: 'json',
    }).done(function (response) {

        const desc = response['description']
        const genes = response['genes']
        
        document.getElementById(id).value = genes;
        if (descid != '') {
            document.getElementById(descid).value = desc;
        }
        geneCount(genes, count_id)
    });
}

function fillSet2(id, descid, count_id) {
    $.ajax({
        url: "getexample2",
        type: "POST",
        data: {},
        dataType: 'json',
    }).done(function (response) {

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
//
// NIDDK RESEARCHER CODE 
//
function example_researcher(tag) {
    console.log(tag.text)
    $('#search_researcher').val(tag.text);
    search_researcher(tag.text);
}

function createResearchCardGridElement(filenames){
    var elem = document.createElement('div');
    elem.className = 'grid-item-researcher';
    elem.setAttribute("data-toggle", "modal")
    elem.setAttribute("data-toggle", "imagecardmodal")
    elem.innerHTML = `<img class="cardimage" onclick="showResearcherCardModal(this)" data-src=${filenames[0]} src=${filenames[1]}>`
    elem.style.position = null
    // console.log(elem)
    return elem
}
function search_researcher(name) {
    $('.researcher-grid').empty();
    //Inlcude this grid sizer for proper layout and sizing of the masonry image layout
    $('.researcher-grid').html('<div class="grid-sizer"></div>');
    $('#researcher-cy-network').empty();
    document.getElementById('researcher-network-title').innerText = '';
    document.getElementById('researcher-cy-network').style.display = 'none';
    document.getElementById('researcher-title').innerText = `Research Summary Information for ${name}`;
    document.getElementById('researcher-data-loading-none').style.display = 'none';
    document.getElementById('researcher-data-loading').style.display = 'block';
    document.getElementById('navbar-toc').style.display = 'block'
    document.getElementById('results_researcher').style.display = 'block'
    //Api call to the backend in order to load the masonry data images for the researcher. Need to implement the API calls if their
    //data is not there.
    const formData = JSON.stringify({'researcher':name})
    console.log(formData)
    fetch("researcher_name", {
        method: "POST",
        headers: {
            'Accept': 'application/json',
            'Content-Type':'application/json'
        },
        body: formData,
    }).then(response => response.json())
    .then(data_dict => {
        const file_names = data_dict['png_paths']
        if (file_names.length == 0) {
            console.log('No images for this researcher. Need to do manual API calls.');
            document.getElementById('researcher-data-loading').style.display = 'none';
            document.getElementById('researcher-data-loading-none').style.display = 'block';
            document.getElementById('researcher-data-loading-none').innerHTML = "NOT IN DATABASE. QUERYING INFORMATION FROM API CALLS. NEED TO IMPLEMENT";
        }else{
            elements = file_names.map(filename => createResearchCardGridElement(filename))
            console.log(elements)
            var $grid = $('.researcher-grid').masonry({
                itemSelector: '.grid-item-researcher',
                columnWidth: '.grid-sizer'
              });
            console.log($grid)
            var $elems = $( elements )
            $grid.append( $elems ).masonry( 'appended', $elems );
            document.getElementById('researcher-data-loading').style.display = 'none'
        }
    })

    //MAKING API CALL TO NIDDK KG and Loading base Cytoscape
    route_for_kg_connections = `http://localhost:3000/api/knowledge_graph?start=Principal Investigator&start_term=${name}&start_field=label&limit=25&end=Principal Investigator&relation=PI Gene,PI Diseases`
    fetch(route_for_kg_connections, {
        method: "GET",
        headers: {
            'Accept': 'application/json',
            'Content-Type':'application/json', 
        }
    }).then(response => response.json())
    .then((data) => {
        console.log(typeof(data))
        if (data.length == 0){
            document.getElementById('researcher-network').innerText = "No connections to other Researchers";
        }else{
            
            document.getElementById('researcher-network-title').innerText = `Researcher Subnetwork for ${name}`;
            text_content = ``
            for (i = 0; i< data.length; i++) {
                console.log(data[i]['data'])
                text_content += `<p>${JSON.stringify(data[i]['data'])}</p>`
            }
            document.getElementById('researcher-network').innerText = text_content;

            document.getElementById('researcher-cy-network').style.display = 'block';
            network_div = document.getElementById('researcher-cy-network');
            var cy = cytoscape({

                container: network_div, // container to render in
              
                elements: data,
              
                style : [{
                    "selector": 'node',
                    "style": {
                    'background-color': 'data(color)',
                    'border-color': 'data(borderColor)',
                    'border-width': 'data(borderWidth)',
                    'label': 'data(label)',
                    "text-valign": "center",
                    "text-halign": "center",
                    'width': "50",
                    'height': "50",
                    }
                },
                {
                    "selector": 'edge',
                    "style": {
                    'curve-style': 'straight',
                    'line-color': 'data(lineColor)',
                    'width': '3',
                    'label': 'data(relation)',
                    "text-rotation": "autorotate",
                    "text-margin-x": "0px",
                    "text-margin-y": "0px",
                    'font-size': '12px',
                    'target-arrow-shape': "data(directed)",
                    'target-endpoint': 'outside-to-node',
                    'source-endpoint': 'outside-to-node',
                    'target-arrow-color': 'data(lineColor)',
                    }
                },
                {
                    "selector": 'node.highlight',
                    "style": {
                        'border-color': 'gray',
                        'border-width': '2px',
                        'font-weight': 'bold',
                        'font-size': '18px',
                        'width': "90",
                        'height': "90",
                    }
                },
                {
                    "selector": 'node.focused',
                    "style": {
                        'border-color': 'gray',
                        'border-width': '2px',
                        'font-weight': 'bold',
                        'font-size': '18px',
                        'width': "90",
                        'height': "90",
                    }
                },
                {
                    "selector": 'edge.focusedColored',
                    "style": {
                        'line-color': '#F8333C',
                        'width': '6'
                    }
                },
                {
                    "selector": 'node.semitransp',
                    "style":{ 'opacity': '0.5' }
                },
                {
                    "selector": 'node.focusedSemitransp',
                    "style":{ 'opacity': '0.5' }
                },
                {
                    "selector": 'edge.colored',
                    "style": {
                        'line-color': '#F8333C',
                        'target-arrow-color': '#F8333C',
                        'width': '6'
                    }
                },
                {
                    "selector": 'edge.semitransp',
                    "style":{ 'opacity': '0.5' }
                },
                {
                    "selector": 'edge.focusedSemitransp',
                    "style":{ 'opacity': '0.5' }
                }],
              
                layout: {
                  name: 'breadthfirst'
                }
              
              });

        }

    })

}
function showResearcherCardModal(element) {
    // The element parameter here is the img tag with the src holding the path to the image. 
    image_path = element.getAttribute("data-src")
    // console.log(image_path)
    modal_body = document.getElementById('modalbodycontent')
    // console.log(modal_body)
    modal_body.innerHTML = ''
    modal_body.innerHTML = `<img src="${image_path}" style=height:400px;>`
    $('#imagecardmodal').modal()
}

// Autocomplete for the input for NIDDK researcher page
const researcher_list = fetch("static/researcher_list.json").then(data => {
    return data.json()

}).then(data => data);

(function () {
    $('#search_researcher').val('');
    const autoCompleteJSResearcher = new autoComplete({
        selector: "#search_researcher",
        placeHolder: "Input an NIDDK researcher name",
        data: {
            src: researcher_list
        },
        resultsList: {
            element: (list, data) => {
                const info = document.createElement("p");
                if (data.results.length > 0) {
                    info.innerHTML = `Displaying <strong>${data.results.length}</strong> out of <strong>${data.matches.length}</strong> results`;
                } else {
                    info.innerHTML = `Found <strong>${data.matches.length}</strong> matching results for <strong>"${data.query}"</strong>`;
                }
                list.prepend(info);
            },
            noResults: true,
            maxResults: 15,
            tabSelect: true
        },
        resultItem: {
            highlight: true
        },
        events: {
            input: {
                selection: (event) => {
                    autoCompleteJSResearcher.input.value = event.detail.selection.value;
                    search_researcher(event.detail.selection.value);
                }
            }
        }
    });
})();


// //Map implementation using leaflet and tile maps from online
// var map = L.map('map').setView(['37.0902','-95.7129'], 4);


// // L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
// // maxZoom: 7,
// // attribution: '© OpenStreetMap'
// // }).addTo(map);

// L.tileLayer('https://{s}.tile.jawg.io/jawg-sunny/{z}/{x}/{y}{r}.png?access-token={accessToken}', {
// 	attribution: '<a href="http://jawg.io" title="Tiles Courtesy of Jawg Maps" target="_blank">&copy; <b>Jawg</b>Maps</a> &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
// 	minZoom: 0,
// 	maxZoom: 15,
// 	subdomains: 'abcd',
// 	accessToken: 'kP5zFalUlpUkkreJz6QYMRiXibhtyZVMNYoOv199RUSjrJf73RCdSjXYRxjucPgk'
// }).addTo(map);


// async function getLocData(){
//     const location_info = await fetch("static/location_info.json").then(data => data.json()).then(data => data);
//     let count_here = 0
//     for (const key in location_info) {
//         console.log('Here')
//         count_here += 1

//         html_string = ''
//         html_string = '<div class="scroll">'
//         for (i= 0 ; i < location_info[key].length; i++) {
//             var marker = L.marker(key.split(',')).addTo(map);
//             html_string += `<a href=#search_researcher onClick="example_researcher(this)">${location_info[key][i]}</a><br>`
//         }
//         html_string += '</div>'
//         marker.bindPopup(html_string);
//         if (count_here == 30){
//             break
//         }

//     }
//     console.log('In function get loc')

// }

// getLocData()

function toggle_map(){
    map_div = document.getElementById('map-div-d3');
    map_button = document.getElementById('map-button');
    if (map_div.style.display == 'none') {
        map_div.style.display = 'block';
    }
    else{
        map_div.style.display = 'none';
    }
    console.log(map_button.innerText);
    if (map_button.innerText.includes('Hide')) {
        map_button.innerText = "Show Map";
    }
    else{
        map_button.innerText = "Hide Map";
    }
}


// async function makeLocMapSVG(){
//     const location_info = await fetch("static/location_info_d3.json").then(data => data.json()).then(data => data);
//     const margin = { top: 5, right: 5, bottom: 5, left: 5 };
// 	let width = document.querySelector("body").clientWidth;
// 	let height = 700;
//     let projection = d3.geoEquirectangular().center([0, 0]);
//     const pathGenerator = d3.geoPath().projection(projection);

//     const hexbin = d3.hexbin()
//     .extent([[0, 0], [width, height]])
//     .radius(3)
//     .x(d => d.xy[0])
//     .y(d => d.xy[1]);

//     const bins = hexbin(location_info.map(d => ({xy: projection([d.longitude, d.latitude])})))

//     const svg = d3.select("#map-svg").attr("viewBox", [0, 0, width, height]);

//     let g = svg.append("g");
//     d3.json(
//             "https://raw.githubusercontent.com/iamspruce/intro-d3/main/data/countries-110m.geojson"
//         )
//         .then((data) => {
//             g.selectAll("path").data(data.features).join("path").attr("d", pathGenerator);
//             svg.append('g').selectAll("path")
//             .data(bins)
//             .join("path")
//             .attr("transform", d => `translate(${d.x},${d.y})`)
//             .attr("d", d => hexbin.hexagon())
//             .attr("fill", 'green')
//             .attr("stroke", 'white');
//         });

// }

// async function makeUSALocMapSVG(){
//     const location_info = await fetch("static/location_info_d3.json").then(data => data.json()).then(data => data);
//     const margin = { top: 5, right: 5, bottom: 5, left: 5 },
// 	width = document.querySelector("body").clientWidth,
// 	height = 700;
//     let projection = d3.geoAlbers();
//     const pathGenerator = d3.geoPath().projection(projection);

//     const hexbin = d3.hexbin()
//     .extent([[0, 0], [width, height]])
//     .radius(10)
//     .x(d => d.xy[0])
//     .y(d => d.xy[1]);
//     const bins = hexbin(location_info.map(d => ({xy: projection([d.longitude, d.latitude])})))

//     const svg = d3.select("#map-svg2").attr("viewBox", [0, 0, width, height]);
//     let g = svg.append("g");
//     d3.json("static/counties-10m.json")
//     .then((data) => {
//             stateMesh = topojson.mesh(data, data.objects.states)
//             g.append("path")
//             .datum(stateMesh)
//             .attr("fill", "none")
//             .attr("stroke", "#777")
//             .attr("stroke-width", 0.5)
//             .attr("stroke-linejoin", "round")
//             .attr("d", d3.geoPath(projection));
//             // g.append('path').attr('fill', )

//             svg.append('g').selectAll("path")
//             .data(bins)
//             .join("path")
//             .attr("transform", d => `translate(${d.x},${d.y})`)
//             .attr("d", d => hexbin.hexagon())
//             .attr("fill", 'black')
//             .attr("stroke", 'white');
//         });

// }

async function makeUSAMap(){
    document.getElementById("map-loading").innerHTML = "<div class='loader justify-content-center'></div>";
    console.log('Width of site body')
    console.log(document.querySelector("body").clientWidth);
    const width = 975;
    const height = 610;
    //Setting up the projection in order to convert coordinates into points that fit on the map in screen
    const projection = d3.geoAlbers().scale(1300).translate([487.5,305]);
    const path  = d3.geoPath();
    const svg = d3.select("#map-svg3").attr("height", height).attr('width', '100%').attr("viewBox", [0, 0, width, height]);
    //.style('max-width', '100%').style('height', 'auto');

    const us = await fetch("static/states-albers-10m.json").then(data => data.json()).then(data => data);
    //Creating the background and outline of the US using topojson module and the path attribute for svg elements
    const statesBackground = svg
    .append('path')
    .attr('fill', '#ddd')
    .attr('d', path(topojson.feature(us, us.objects.nation)));
    //Appnding another path element to make the shape of each of the states
    const statesBorders = svg
    .append('path')
    .attr('fill', 'none')
    .attr('stroke', '#fff')
    .attr('stroke-linejoin', 'round')
    .attr('stroke-linecap', 'round')
    .attr('d', path(topojson.mesh(us, us.objects.states, (a, b) => a !== b)));

    let location_info = await fetch("static/location_info_d3.json").then(data => data.json()).then(data => data);
    location_info = location_info.slice(0, 3);
    console.log(location_info)
    const researcherLocationElements = svg
    .selectAll('g')
    .data(location_info)
    .join('g');

    const researcherLocationPoints = researcherLocationElements
    .append('g')
    .attr('class', 'researcherLocationG')
    .attr('transform', (d) => `translate(${projection([d.longitude, d.latitude]).join(",")})`);
    
    researcherLocationPoints.append('circle')
    .attr('r', 5)
    .attr('class', 'researcherLocationPoint')
    .attr('stroke', 'white')
    .style('fill', 'black');


    // let tooltip = d3
    // .select("body")
    // .append("div")
    // .attr("class", "tooltip scroll")
    // .attr("id", "researchertooltiplist")
    // .style('border-style', 'solid')
    // .style('background-color', '#ddd')
    // .style("opacity", 0);

    // d3.selectAll('circle')
    // .on("click", (event, data) => {
    //     tooltip.transition().duration(200).style("opacity", 0.9);
    //     tooltip.style("display", 'block');
    //     html_string = '<button type="button" class="close sticky-top" onClick="closeResearcherList()" ><span aria-hidden="true">&times;</span></button><br>'
    //     for (i= 0 ; i < data.names.length; i++) {
    //         html_string += `<a href=#search_researcher onClick="example_researcher(this)">${data.names[i]}</a><br>`
    //     }
    //     tooltip.html(html_string)
    //     .style("left", event.pageX + "px")
    //     .style("top", event.pageY - 28 + "px");

    // })


    let tooltip = d3
    .select("body")
    .append("div")
    .attr("class", 'tooltip')
    .attr("id", "researchertooltiplist")
    .style('border-style', 'solid')
    .style('border-radius', '10px')
    .style('box-shadow', '3px 4px 5px black')
    .style('background-color', '#ddd')
    .style("opacity", 0);
    d3.selectAll('circle')
    .on("click", (event, data) => {
        html_string = `<button type="button" class="close sticky-top" onClick="closeResearcherList()" ><span aria-hidden="true">&times;</span></button><br>`
        html_string += `<p class=text-center><strong>${data.institution}</strong></p>`
        html_string += `<div class="scroll">`
        document.getElementById('researchertooltiplist').innerHTML = '';
        console.log(data)
        for (i= 0 ; i < data.names.length; i++) {
            html_string += `<a href=#search_researcher onClick="example_researcher(this)">${data.names[i]}</a><br>`
        }
        html_string += `</div>`
        tooltip.html(html_string)
        .style("left", event.pageX + "px")
        .style("top", event.pageY - 28 + "px");
        tooltip.transition().duration(200).style("opacity", 0.9);
        tooltip.style("display", 'block');

    })

    let zooming = d3
    .zoom()
    .scaleExtent([1, 3])
    .on("zoom", (event) => {
        svg.selectAll("path").attr("transform", event.transform);
        console.log(event.transform)
        console.log(event.transform.k)
        console.log(event.transform.x)
        console.log(event.transform.y)
        svg.selectAll("circle")
          .attr("transform", event.transform)
          .attr("r", 6 / event.transform.k);
          
  });
  svg.call(zooming)
  document.getElementById("map-loading").innerHTML = "";

}



// makeUSAMap();


async function makeUSAMapZoom(){
    document.getElementById("map-loading").innerHTML = "<div class='loader justify-content-center'></div>";
    const width = 975;
    const height = 610;
    //Setting up the projection in order to convert coordinates into points that fit on the map in screen
    const projection = d3.geoAlbers().scale(1300).translate([487.5,305]);
    const path  = d3.geoPath();
    const svg = d3.select("#map-svg3").attr("height", height).attr('width', '100%').attr("viewBox", [0, 0, width, height]);
    //.style('max-width', '100%').style('height', 'auto');

    const us = await fetch("static/states-albers-10m.json").then(data => data.json()).then(data => data);
    //Creating the background and outline of the US using topojson module and the path attribute for svg elements
    var g = svg.append("g");
    const statesBackground = g
    .append('path')
    .attr('fill', '#ddd')
    .attr('d', path(topojson.feature(us, us.objects.nation)));

    //Appnding another path element to make the shape of each of the states
    const statesBorders = g
    .append('path')
    .attr('fill', 'none')
    .attr('stroke', '#fff')
    .attr('stroke-linejoin', 'round')
    .attr('stroke-linecap', 'round')
    .attr('d', path(topojson.mesh(us, us.objects.states, (a, b) => a !== b)));

    let location_info = await fetch("static/location_info_d3.json").then(data => data.json()).then(data => data);
    // location_info = location_info.slice(0, 3);
    console.log(location_info);


    
    var circle = svg
    .select("g").selectAll("circle")
    .data(location_info);


    circle.enter().append("circle")
    .attr("class", "researcherLocationPoint")
    .attr("cx", function(d){ return projection([d.longitude, d.latitude])[0] })
    .attr("cy", function(d){ return projection([d.longitude, d.latitude])[1] })
    .attr("r", 4)
    .style("fill", "black")
    .attr("stroke", "#012B4E")
    .style("opacity", 0.9);

    let tooltip = d3
    .select("body")
    .append("div")
    .attr("class", 'tooltip')
    .attr("id", "researchertooltiplist")
    .style('border-style', 'solid')
    .style('border-radius', '10px')
    .style('box-shadow', '3px 4px 5px black')
    .style('background-color', '#ddd')
    .style("opacity", 0);

    d3.selectAll('circle')
    .on("click", (event, data) => {
        html_string = `<button type="button" class="close sticky-top" onClick="closeResearcherList()" ><span aria-hidden="true">&times;</span></button><br>`
        html_string += `<p class=text-center><strong>${data.institution}</strong></p>`
        html_string += `<div class="scroll">`
        document.getElementById('researchertooltiplist').innerHTML = '';
        console.log(data)
        for (i= 0 ; i < data.names.length; i++) {
            html_string += `<a href=#search_researcher onClick="example_researcher(this)">${data.names[i]}</a><br>`
        }
        html_string += `</div>`

        tooltip.html(html_string)
        .style("left", event.pageX + "px")
        .style("top", event.pageY - 28 + "px");

        tooltip.transition().duration(200).style("opacity", 0.9);
        tooltip.style("display", 'block');

    })

    let zooming = d3
    .zoom()
    .scaleExtent([1, 50])
    .on("zoom", (event) => {
        g.attr("transform", event.transform);
        console.log(event.transform)
        svg.selectAll("circle")
          .attr("r", 4 / event.transform.k);
          
  });

  function recenter() {
    console.log('recentering')
    try {
      // get a handle to the transform object, and reset it to the
      // zoom identity (resets to zoom out all the way)
      let d3zoom = zooming
      let ztrans = d3zoom.transform
      let t = d3.zoomIdentity 
      svg.transition().duration(1000).call(ztrans,t)

      // get the object and measurements to determine the virtual
      // "bounding" rectangle around the directional extremes of sites
      // i.e., north (topmost), east (rightmost), etc
      let circles    = svg.select('g').selectAll('circle')
      let data       = circles.data()
      let long       = data.map(o => parseFloat(o.longitude))
      let lat        = data.map(o => parseFloat(o.latitude))
      let rightmost  = d3.max(long)
      let leftmost   = d3.min(long)
      let topmost    = d3.max(lat)
      let bottommost = d3.min(lat)
      // convert the lat/long to points in the projection
      let lt = projection([leftmost,topmost])
      let rb = projection([rightmost,bottommost])
      // calc the gaps (east - west, south - north) in pixels
      let g  = rb[0]-lt[0]
      let gh = rb[1]-lt[1]

      // get the dimensions of the panel in which the map sits
      let w  = svg.node().parentElement.getBoundingClientRect().width
      let h  = svg.node().parentElement.getBoundingClientRect().height

      // the goal here is to move the halfway point between leftmost and rightmost
      // on the projection sites to the halfway point of the panel in which the
      // svg element resides, and same for 'y'
      // the scale is the value 90% of the scale factor of either east-west to width
      // or south-north to height, whichever is greater
      let neoScale = 0.9 / Math.max(g/w, gh/h)

      // now recalculate what will be the difference between the current
      // center and the center of the scaled virtual rectangle
      // this finds the difference between the centers
      // the new center of the scaled rectangle is the average of the left and right
      // or top and bottom points
      let neoX = w/2 - (neoScale * ((lt[0]+rb[0])/2))
      let neoY = h/2 - (neoScale * ((lt[1]+rb[1])/2))

      // TRANSLATE FIRST!  then scale.
      t = d3.zoomIdentity.translate(neoX,neoY).scale(neoScale)
      svg.transition().duration(1000).call(ztrans, t)
    }
    catch(e)
    {
      console.log(e)
    }
  }


svg.call(zooming);


document.getElementById("map-loading").innerHTML = "";



}

makeUSAMapZoom();

function closeResearcherList(){
    document.getElementById('researchertooltiplist').innerHTML = ''
    document.getElementById('researchertooltiplist').style.display = 'none';
}
$(document).ready(function () {

    // SMALL NAV MENU

    $("<select class='selectize-nav text-center justify-content-center' id='hamburger'> <img href='static/img/hamburger.jpeg/> </select>").appendTo("#mainnav");

    // Create default option "Go to..."
    $("<option />", {
        "selected": "selected",
        "value": "",
        "text": "Go to.."
    }).appendTo("nav select");

    // Populate dropdown with menu items
    $("nav a").each(function () {
        var el = $(this);
        if (el.attr("href").substr(0, 1) != "#" && el.attr("href") != '/' && el.attr("id") != 'toc') {
            $("<option />", {
                "value": el.attr("href"),
                "text": el.text().trim()
            }).appendTo("#mainnav select");
        }
    });

    $("nav select").change(function () {
        window.location = $(this).find("option:selected").val();
    });

    // BOLD CURRENT PAGE

    $('.nav-link').each(function () {
        var url = window.location.href
        if (url.includes('#')) {
            url = url.split('#')[0]
        }

        if ($(this).prop('href') == url) {
            $(this).addClass('active');
            $(this).parents('li').addClass('active');
        }
        // if (url.split('/')[3].startsWith('GSE')) {
        //     $("#viewer").addClass('active'); 
        //     $("#viewer").parents('li').addClass('active');
        // }
    });

    var currURL = window.location.href.split("/");


    var first_load = true;

    // SELECTIZE FOR VIEWER AND HOME PAGE APPYTER T2D

    var $gene_select = $('#gene-select').selectize({
        preload: true,
        valueField: 'gene_symbol',
        labelField: 'gene_symbol',
        searchField: 'gene_symbol',
        render: {
            option: function (item, escape) {
                return '<div class="pt-2 light">' + item.gene_symbol + '</div>';
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
                        first_load = false;
                    } else if (first_load) {
                        $gene_select[0].selectize.setValue(res[0]['gene_symbol']);
                        first_load = false;
                    }
                }
            });
        },
        persist: true,
    });


    // Configure dropdown menu

    $('.dropdown-menu a.dropdown-toggle').on('click', function (e) {
        if (!$(this).next().hasClass('show')) {
            $(this).parents('.dropdown-menu').first().find('.show').removeClass('show');
        }
        var $subMenu = $(this).next('.dropdown-menu');
        $subMenu.toggleClass('show');


        $(this).parents('li.nav-item.dropdown.show').on('hidden.bs.dropdown', function (e) {
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
        $('#boxplotloader').addClass('loader');

        // Gene
        var gene_symbol = $('#gene-select').val();

        // Conditions FOR SINGLE CELL NEED TO CHANGE BUT FOR CLUSTERS NOW WHAT ARE THE SELECTED CLUSTERS
        var conditions = [];
        $('.condition-btn.plotted').each(function () { conditions.push($(this).attr('data-group_label')) }); conditions

        // AJAX Query
        $.ajax({
            url: $('#boxplot').attr('data-url-plot'), //"{{ url_for('plot_api') }} + "/" + $('#boxplot').attr('data-geo-acc'),
            method: 'post',
            data: JSON.stringify({ 'gene_symbol': gene_symbol, 'conditions': conditions }),
            contentType: 'application/json',
            dataType: 'json',
            error: function () {
                "<p> Gene not found <\p>"
            },
            success: function (res) {
                $('#boxplotloader').removeClass('loader');
                layout = {
                    plot_bgcolor: "#00FFFFFF",
                    paper_bgcolor: "#00FFFFFF"
                }
                Plotly.newPlot('boxplot', res['data'], res['layout'], config = { responsive: true }); // maybe plotly.react will be faster here
            }
        });

    };


    // 2. Listeners
    // Gene
    var is_gse = currURL.filter(x => x.includes('GSE'))
    if (is_gse.length > 0) {
        var boxplot_selectize = $gene_select[0].selectize;
        boxplot_selectize.on('change', function (value) {
            boxplot();
        })
    }

    // Conditions
    $('.condition-btn').on('click', function (evt) {
        $(this).toggleClass('plotted'); // making a specific button plotted or not
        boxplot();
    })

    // MAKE STUIDIES TABLE A DataTable

    $('#studies-table').DataTable({
        dom: 'Bfrtip',
        buttons: [
            'copy', { extend: 'csv', title: `D2H2-studies-table` }
        ],
        responsive: true
    });


    // OPEN CUSTOMIZED WORKFLOW DEPENDING ON SELECTION IN SCG

    $('#scg-link').on('click', function () {

        var workflow = document.getElementById("workflow").value;


        if (workflow) {

            $('#scg-link').prop('href', workflow);
            return;
        }
        $('#scg-link').prop('href', "https://scg.maayanlab.cloud/");
    });

    // OPEN TO DEG IN BULK RNA SEQ ANALYSIS
    $('#dgea-button').on('click', async function () {
        var control_condition = $('#condition-select-control').val();
        var perturb_condition = $('#condition-select-perturb').val();
        var species = document.getElementById("species").innerText

        if (!control_condition || !perturb_condition) {
            alert("Please select both conditions")
            return;
        }

        if (control_condition == perturb_condition) {
            alert("Please select two different conditions")
            return;
        }

        document.getElementById("dgea-loading").innerHTML = "<div class='loader justify-content-center'></div>";

        var gse = document.getElementById("gse").innerText

        var gsedata = JSON.stringify({ 'gse': gse, 'control': control_condition, 'perturb': perturb_condition, 'species': species });


        $.ajax({
            url: "api/data",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done(async function (response) {

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

            var download_link = "https://appyters.maayanlab.cloud/Bulk_RNA_seq/" + id.session_id + "/DEG_results_" + control_condition.split(" ").join('%20') + "%20vs.%20" + perturb_condition.split(" ").join('%20') + ".csv"


            const final_url = "https://appyters.maayanlab.cloud/Bulk_RNA_seq/" + id.session_id + "/#differential-gene-expression"
            window.open(final_url, '_blank')

        });
    });

    // PERFORM DIFFERENTIAL GENE ANALYSIS AND CREATE RELEVANT TABLE
    $('#dge-button').on('click', async function () {
        clear_dge()
        var control_condition = $('#condition-select-control').val();
        var perturb_condition = $('#condition-select-perturb').val();

        if (!control_condition || !perturb_condition) {
            alert("Please select both conditions")
            return;
        }

        if (control_condition == perturb_condition) {
            alert("Please select two different conditions")
            return;
        }

        document.getElementById("dge-loading").innerHTML = "<div class='loader justify-content-center'></div> <div class='text-center'><p> It may take 2-3 minutes to compute DEGs</p><div>";

        var gse = document.getElementById("gse").innerText
        var species = document.getElementById("species").innerText
        var method = document.getElementById("method").value

        var logCPM = document.getElementById("logCPM").checked
        var log = document.getElementById("log").checked
        var q = document.getElementById("q").checked
        var z = document.getElementById("z").checked
        var norms = { 'logCPM': logCPM, 'log': log, 'q': q, 'z': z }


        var gsedata = JSON.stringify({ 'gse': gse, 'control': control_condition, 'perturb': perturb_condition, 'method': method, 'species': species, 'norms': norms });
        $.ajax({
            url: "dgeapi",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done(async function (response) {
            document.getElementById("dge-loading").innerHTML = "";
            var plot = response['plot']
            var table = response['table']

            var rows = table.split('\n').slice(1, -1);
            clear_dge()
            window.Bokeh.embed.embed_item(plot)

            var table_id = 'dge-table'
            var tabletext = `<table id='${table_id}' class='styled-table'><thead><tr>`
            const clear_button = "<div class='mx-auto justify-content-center text-center'><button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_dge();'> Clear Results </button></div>"
            var url = currURL.filter(x => !x.includes('GSE')).join('/') + '/singlegene'
            var pvals;
            var adjpvals;
            var logfc;

            if (method === 'limma') {
                tabletext += "<th></th><th>Adj. P Value</th><th>P Value</th><th>t</th><th>AvgExpr</th><th>logFC</th><th>B</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[5]).toPrecision(4) + "</td><td>" + Number(vals[4]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td><td>" + Number(vals[1]).toPrecision(4) + "</td><td>" + Number(vals[6]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
                logfc = table.column(5).data()

            } else if (method === 'characteristic_direction') {
                tabletext += "<th></th><th>CD Coeffecient</th><th>P-value</th><th>LogFC</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[1]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(2).data()
                pvals = table.column(2).data()
                logfc = table.column(3).data()

            } else if (method === 'DESeq2') {
                tabletext += "<th></th><th>Adj. P-value</th><th>P-value</th><th>lfcSE</th><th>stat</th><th>baseMean</th><th>log2FC</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[6]).toPrecision(4) + "</td><td>" + Number(vals[5]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td><td>" + Number(vals[4]).toPrecision(4) + "</td><td>" + Number(vals[1]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
                logfc = table.column(6).data()
            }
            var genes = table.column(0).data()


            genes = genes.map(x => x.replace(/<\/?[^>]+(>|$)/g, ""))

            var genelist_buttons =
            `<div class="row justify-content-center mx-auto text-center">
                <div class="h7">Submit the top</div>
                <input class="" id='numgenes' type='number' step='1' value='100' pattern='[0-9]' min='1' class='m-2'
                    style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;' />
                <select id='dir' style='margin-left: 10px; margin-right: 10px;'>
                    <option value='up'>upregulated</option>
                    <option value='down'>downregulated</option>
                </select>
                <div class="h7"> differentially expressed genes with a
                    <select id='col-to-use' style='margin-left: 10px; margin-right: 10px;'>
                        <option value='pval'>p-value</option>
                        <option value='adjpval'>adjusted p-value</option>
                    </select>
                </div>
                <div class=" h7"> less than </div>
                <input class='' id='signifigance' type='number' step='.01' value='.05' max='1'
                    style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;"' />
                <div class="h7"> to</div>
            </div>
            <div class="row justify-content-center mx-auto text-center">
                <button class="btn btn-primary btn-group-sm m-2"
                    onclick="submit_geneset('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${logfc.join(',')}')">Gene Set Queries</button>
            </div>
            <div class="row justify-content-center mx-auto text-center">
                <button type="button" class="btn btn-primary btn-group-sm m-2"
                    onclick="filter_and_submit_to_enrichr('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${logfc.join(',')}', '${gse}-${control_condition}-vs-${perturb_condition}')">
                    Enrichr
                    <img src="/static/img/enrichrlogo.png" class="img-fluid mr-3" style="width: 45px" alt="Enrichr">
                </button>
                <button type="button" class="btn btn-primary btn-group-sm m-2"
                    onclick="filter_and_submit_to_kg('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${logfc.join(',')}', '${gse}-${control_condition}-vs-${perturb_condition}')">
                    Enrichr-KG
                    <img src="/static/img/enrichr-kg.png" class="img-fluid mr-3" style="width: 45px" alt="Enrichr">
                </button>
            </div>
            `

            document.getElementById("geneset-buttons").innerHTML = (clear_button + genelist_buttons)
        })
            .fail(function (jqXHR, textStatus, errorThrown) {
                alert('An internal server error occured, please try again')
                document.getElementById("dge-loading").innerHTML = "";
            })
    });



    // PERFORM DIFFERENTIAL GENE ANALYSIS AND CREATE RELEVANT TABLE
    $('#dge-button-single').on('click', async function () {
        clear_dge_single()

        document.getElementById("dge-loading").innerHTML = "<div class='loader justify-content-center'></div>";

        var gse = document.getElementById("singlegse").innerText
        var species = document.getElementById("species").innerText
        var method = document.getElementById("method").value
        var condition_group = document.getElementById("methodsingle").value
        var diffcluster = document.getElementById("differentialcluster").value

        var logCPM = document.getElementById("logCPM").checked
        var log = document.getElementById("log").checked
        var q = document.getElementById("q").checked
        var z = document.getElementById("z").checked
        var norms = { 'logCPM': logCPM, 'log': log, 'q': q, 'z': z }

        var gsedata = JSON.stringify({ 'gse': gse, 'species': species, 'conditiongroup': condition_group, 'method': method, 'norms': norms, 'diffcluster': diffcluster });

        $.ajax({
            url: "/dgeapisingle",
            contentType: 'application/json',
            type: "POST",
            dataType: 'json',
            data: gsedata
        }).done(async function (response) {
            document.getElementById("dge-loading").innerHTML = "";
            var plot = response['plot']
            var table = response['table']
            var desc = response['description']

            var rows = table.split('\n').slice(1, -1);
            clear_dge_single()
            window.Bokeh.embed.embed_item(plot)

            var table_id = 'dge-table'
            var tabletext = `<table id='${table_id}' class='styled-table'><thead><tr>`
            const clear_button = "<div class='mx-auto justify-content-center text-center'><button type='button' class='btn btn-dark btn-group-sm mt-3 mb-3' onclick='clear_dge_single();'> Clear Results </button></div>"
            var url = currURL.filter(x => !x.includes('GSE')).join('/') + '/singlegene'
            var pvals;
            var adjpvals;
            var logfc;

            if (method === 'limma') {
                tabletext += "<th></th><th>Adj. P Value</th><th>P Value</th><th>t</th><th>AvgExpr</th><th>logFC</th><th>B</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[5]).toPrecision(4) + "</td><td>" + Number(vals[4]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td><td>" + Number(vals[1]).toPrecision(4) + "</td><td>" + Number(vals[6]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
                logfc = table.column(5).data()

            } else if (method === 'characteristic_direction') {
                tabletext += "<th></th><th>CD Coeffecient</th><th>P-value</th><th>LogFC</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[1]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(2).data()
                pvals = table.column(2).data()
                logfc = table.column(3).data()

            } else if (method === 'DESeq2') {
                tabletext += "<th></th><th>Adj. P-value</th><th>P-value</th><th>lfcSE</th><th>stat</th><th>baseMean</th><th>log2FC</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[6]).toPrecision(4) + "</td><td>" + Number(vals[5]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td><td>" + Number(vals[4]).toPrecision(4) + "</td><td>" + Number(vals[1]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
                logfc = table.column(6).data()

            } else if (method === 'wilcoxon') {
                tabletext += "<th></th><th>Adj. P-value</th><th>P-value</th><th>logfoldchanges</th><th>scores</th></tr><tbody>"
                rows.forEach(function (row) {
                    var vals = row.replace(/\s\s+/g, ' ').split(' ');
                    tabletext += `<tr><td><a onclick="setGene('${vals[0]}')" href='${url}' target='_blank'>` + vals[0] + "</a></td><td>" + Number(vals[4]).toPrecision(4) + "</td><td>" + Number(vals[3]).toPrecision(4) + "</td><td>" + Number(vals[2]).toPrecision(4) + "</td><td>" + Number(vals[1]).toPrecision(4) + "</td></tr>"
                });
                tabletext += "</tbody></table>";

                document.getElementById("dge-table-area").innerHTML += tabletext
                var table = $(`#${table_id}`).DataTable({
                    order: [[1, 'asc']],
                    dom: 'Bfrtip',
                    buttons: [
                        'copy', { extend: 'csv', title: name }
                    ]
                });
                adjpvals = table.column(1).data()
                pvals = table.column(2).data()
                logfc = table.column(3).data()
            }
            var genes = table.column(0).data()


            genes = genes.map(x => x.replace(/<\/?[^>]+(>|$)/g, ""))


            document.getElementById("geneset-buttons").innerHTML = (clear_button)
            document.getElementById("enrichment-area").innerHTML = `
            <div class="h4 pl-2 mt-4 mb-4 text-center">Enrichment Analysis for Highly Expressed Genes in ${desc}</div>
            <div class="row justify-content-center mx-auto text-center">
                <div class="h7">Submit the top</div>
                <input class="" id='numgenes' type='number' step='1' value='100' pattern='[0-9]' min='1' class='m-2'
                    style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;' />
                <select id='dir' style='margin-left: 10px; margin-right: 10px;'>
                    <option value='up'>upregulated</option>
                    <option value='down'>downregulated</option>
                </select>

                <div class="h7"> differentially expressed genes with a </div>
                <select id='col-to-use' style='margin-left: 10px; margin-right: 10px;'>
                    <option value='pval'>p-value</option>
                    <option value='adjpval'>adjusted p-value</option>
                </select>
                <div class=" h7"> less than </div>
                <input class='' id='signifigance' type='number' step='.01' value='.05' max='1'
                    style='width: 50px; height: 30px; margin-left: 10px; margin-right: 10px;"' />
                <div class="h7"> to</div>
            </div>
            
            <div class="row justify-content-center mx-auto text-center">
                <button type="button" class="btn btn-primary btn-group-sm mt-3 mb-3"
                    onclick="filter_and_submit_to_enrichr('${genes.join(',')}', '${adjpvals.join(',')}', '${pvals.join(',')}', '${logfc.join(',')}', 'Highly Expressed Genes in ${desc}')">
                    Submit to Enricher
                    <img src="/static/img/enrichrlogo.png" class="img-fluid mr-3" style="width: 45px" alt="Enrichr">
                </button>
            </div>
            `
        })
            .fail(function (jqXHR, textStatus, errorThrown) {
                alert('An internal server error occured, please try again')
                document.getElementById("dge-loading").innerHTML = "";
            })
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
                return '<div class="pt-2 light">' + item.Condition + '</div>';
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
                return '<div class="pt-2 light">' + item.Condition + '</div>';
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

    $('#species-val').on('change', function () {
        if ($(this).is(':checked')) {
            switchStatus = $(this).is(':checked');
            $gene_select[0].selectize.clearOptions();
            $gene_select[0].selectize.load(function (callback) {
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
            $gene_select[0].selectize.load(function (callback) {
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


    //This listener is for when we click on a different group of sample individuals to look at for a study within the single cell viewer on the dropdown select menu. 
    $('#methodsingle').on('change', async function () {
        clear_dge_single()
        document.getElementById("change-loading").innerHTML = "<div class='loader justify-content-center'></div>";
        var gse = document.getElementById("singlegse").innerText
        var species = document.getElementById("species").innerText
        var condition_group = document.getElementById("methodsingle").value
        var gsedata = JSON.stringify({ 'gse': gse, 'species': species, 'conditiongroup': condition_group });
        generate_single_plots()

        $.ajax({
            url: "/getclusterdata",
            contentType: 'application/json',
            type: "POST",
            data: gsedata,
            dataType: 'json'
        }).done(async function (response) {

            const classes = response['classes']
            const metadict = response['metadict']

            document.getElementById("singlecell").innerHTML = ''
            document.getElementById("differentialcluster").innerHTML = ''
            tabletext = ''
            diffselectiontext = ''
            for (var k = 0; k < classes.length; k++) {
                tabletext += `<tr>
                            <td class="p-0 ml-3 border"><button
                                    style="background-color: #ECEFF1; color: black; width: max-content;"
                                    class="btn m-0 rounded-0 py-0 condition-btn active plotted " data-toggle="button"
                                    autocomplete="off"
                                    data-group_label="${classes[k]}">${classes[k]}</button>
                            </td>
    
                            <td class="border" style="padding: 0px 18px; width: max-content;">
                                <button
                                    style="width: max-content;"
                                    class="btn-custom btn-group-sm btn-collapse collapsed d-flex align-items-center text-center"
                                    data-toggle="collapse" data-target="#samples-${gse}" aria-expanded="false"
                                    aria-controls="samples-${gse}">
                                    <div class="text">${metadict[classes[k]]} Cells
                                    </div>
                                </button>
    
    
                            </td>
                            </tr>`
                diffselectiontext += `<option value="${classes[k]}">${classes[k]}</option>`
            }
            document.getElementById("singlecell").innerHTML = tabletext
            document.getElementById("differentialcluster").innerHTML = diffselectiontext


            //Add the listener for the condition buttons back in this function 
            $('.condition-btn').on('click', function (evt) {
                $(this).toggleClass('plotted'); // making a specific button plotted or not
                boxplot();
            })
            //Fill up the gene selectize with the proper gene list based off each study. 
            $gene_select[0].selectize.clearOptions();
            $gene_select[0].selectize.load(function (callback) {
                $.ajax({
                    url: `/api/singlegenes/${gse}/${condition_group}`,
                    dataType: 'json',
                    error: function () {
                        callback();
                    },
                    success: function (res) {

                        callback(res);
                        $gene_select[0].selectize.setValue(res[0]['gene_symbol']);
                        $("#boxplot").attr("data-url-plot", `/api/plot_single/${gse}/${condition_group}`)
                        boxplot()
                    }
                });
            });
        });
        document.getElementById("change-loading").innerHTML = ''
    })


})


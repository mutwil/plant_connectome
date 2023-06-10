/*
This file contains the original JavaScript contained inside gene.html (i.e., inside those <script> tags).
I've taken them out and dumped them into a JavaScript file for the sake of organization (so that the code
doesn't look too bloated and is more maintainable in the long run).
*/

var searchForm = document.getElementById('node-search-form');
searchForm.addEventListener('submit', function (event) {
    event.preventDefault();
    const searchTerm = event.target[0].value.toLowerCase();
    const cy = document.getElementById('cy')._cyreg.cy;
    const matchingNodes = cy.nodes().filter(function(node) {
        const nodeId = node.id().toLowerCase();
        if (nodeId === searchTerm) {
            return true;
        }
        return nodeId.includes(searchTerm);
    })
    cy.nodes().style('opacity', '0.5');
    matchingNodes.style('opacity', '1');
});

// {{ cytoscape_js_code | safe }} --> meant to be put into the script file this originally belonged in.

function recalculateLayout() {
    const cy = document.getElementById('cy')._cyreg.cy;
    const visibleNodes = cy.elements(':visible');
    const layout = visibleNodes.layout({
		name: 'cose',
		animate: true,
		randomize: true,
		idealEdgeLength: 100,
		nodeOverlap: 20,
		refresh: 20,
		fit: true,
		padding: 0,
		boundingBox: undefined,
		nodeDimensionsIncludeLabels: false,
		excludeFromLayout: cy.nodes(':hidden')
	});
    layout.run();
}
	
function downloadTableAsTSV(tableId, filename) {
    const table = document.getElementById(tableId);
    const rows = table.querySelectorAll('tr');
    let tsvContent = '';

    rows.forEach(row => {
        const rowData = [];
        const cells = row.querySelectorAll('th, td');

        cells.forEach(cell => {
            rowData.push(cell.textContent);
        });

        tsvContent += rowData.join('\t') + '\n';
    });

    const blob = new Blob([tsvContent], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);

    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    link.click();
    URL.revokeObjectURL(url);
}

function downloadAsSVG() {
    // Assuming your Cytoscape instance is named 'cy'
    const cy = document.getElementById('cy')._cyreg.cy;

    // Register the SVG Exporter plugin
    cytoscape.use(cytoscapeSvg);

    // Export the network view as an SVG
    const svgContent = cy.svg({ copyStyles: true, bg: 'white' });

    // Modify the downloaded SVG to have black letters
    const svgDOM = new DOMParser().parseFromString(svgContent, 'image/svg+xml');
    const labels = svgDOM.querySelectorAll('text');
    labels.forEach(label => label.setAttribute('fill', '#000000'));
    const modifiedSvgContent = new XMLSerializer().serializeToString(svgDOM);

    // Create a Blob from the SVG content
    const blob = new Blob([modifiedSvgContent], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);

    // Create a link element, set its href to the Blob URL, and trigger a click event to download the SVG
    const link = document.createElement('a');
    link.href = url;
    link.download = 'network.svg';
    link.click();

    // Revoke the Blob URL
    URL.revokeObjectURL(url);
}

document.getElementById('download-pdf').addEventListener('click', downloadAsSVG);

function sortTable(table, column, asc = true) {
    const directionModifier = asc ? 1 : -1;
    const rows = Array.from(table.tBodies[0].querySelectorAll("tr"));

    const sortedRows = rows.sort((a, b) => {
        const aColText = a.querySelector(`td:nth-child(${column + 1})`).textContent.trim();
        const bColText = b.querySelector(`td:nth-child(${column + 1})`).textContent.trim();

      return aColText.localeCompare(bColText) * directionModifier;
    });

    while (table.tBodies[0].firstChild) {
        table.tBodies[0].removeChild(table.tBodies[0].firstChild);
    }
    table.tBodies[0].append(...sortedRows);

    table.setAttribute('data-sort-direction', asc ? 'asc' : 'desc');
    table.setAttribute('data-sort-column', column);
}


document.querySelectorAll('.sortable thead th').forEach((headerCell) => {
    headerCell.addEventListener('click', () => {
        const table = headerCell.parentElement.parentElement.parentElement;
        const columnIndex = Array.from(headerCell.parentElement.children).indexOf(headerCell);
        const currentDirection = table.getAttribute('data-sort-direction') === 'asc' ? true : false;
        sortTable(table, columnIndex, !currentDirection);
    });
});

function showPaperPopup(pubmedID, source, typa, target) {
    // Create the PubMed API URL for the paper
    var pubmedURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=" + pubmedID + "&rettype=abstract&retmode=text";

    // Use jQuery to load the paper data from the PubMed API
    $.get(pubmedURL, function(data) {
    // Create the HTML for the paper popup
        var sourceWords = source.split(/[,\s&+/\\']+/);
        var targetWords = target.split(/[,\s&+/\\']+/);
        var underlinedData = data;
        for (var i = 0; i < sourceWords.length; i++) {
            var sourceRegex = new RegExp("\\b" + sourceWords[i].replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&') + "\\b", "gi");
            underlinedData = underlinedData.replace(sourceRegex, "<u><span class='red-text'>$&</span></u>");
        }
        for (var j = 0; j < targetWords.length; j++) {
            var targetRegex = new RegExp("\\b" + targetWords[j].replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&') + "\\b", "gi");
            underlinedData = underlinedData.replace(targetRegex, "<u><span class='red-text'>$&</span></u>");
        }
        var popupHTML = "<div class='paper-popup'><h2>" + source + " " + typa + " " + target + "</h2><pre>" + underlinedData + "</pre> Click <a href='https://pubmed.ncbi.nlm.nih.gov/" + pubmedID + "'target='_blank' class='pubmed-hyperlink'>here</a> to view the publication.</div>";

        // Display the popup using the jQuery UI dialog function
        var $popup = $(popupHTML).dialog({
            width: "90%",
            maxWidth: "none",
            minHeight: "none",
            modal: true,
            close: function () {
                $popup.dialog("destroy").remove();
                currentPopup = null; // Reset the currentPopup variable
            },
        });

        // Set the currentPopup variable to the new popup
        currentPopup = $popup;

        // Close the popup when clicking outside the box
        $(document).on("click", function (event) {
            if (currentPopup !== null && !$popup.is(event.target) && $popup.has(event.target).length === 0) {
                currentPopup.dialog("destroy").remove();
                currentPopup = null; // Reset the currentPopup variable
            }
        });
    });
}

// Event listener for clicking on pubmed-link elements
$(document).on('click', '.pubmed-link', function() {
    var pubmedID = $(this).data('pubmed-id');
    var source = $(this).data('source');
    var typa = $(this).data('typa');
    var target = $(this).data('target');

    showPaperPopup(pubmedID, source, typa, target);
});

$(window).on("load", function() {
    const cy = document.getElementById("cy_wrapper");
    cy.style.height = `${window.innerHeight * 0.8}px`;

});
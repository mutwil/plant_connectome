// Get forms in order of gene-form, name-form, and title-form:
const forms = document.querySelectorAll('form');

const titleInput = document.getElementById("gene_id");


// Change placeholder text based which button is being hovered over.
const alias_button_events = () => {
    let alias_button = document.querySelector('input[value = Alias]');
    alias_button.addEventListener('mouseover', () => {
        titleInput.placeholder = "e.g., CESA1 (hit 'Enter' to search)";
    })
    alias_button.addEventListener('mouseout', () => {
        titleInput.placeholder = "e.g., CESA (hit 'Enter' to search)";
    })
}

// Attaches event listeners to the gene form.
const gene_form_listeners = () => {
    let form_buttons = forms[0].querySelectorAll('.button');
    let form_actions = (function() {
        let temp = [];
        form_buttons.forEach(button => {
            temp.push(button.getAttribute('formaction'));
        })
        return temp;
    })();
    
    for (let i = 0 ; i < form_buttons.length ; i++) {
        form_buttons[i].addEventListener('click', () => {
            submitGeneForm(event, forms[0], form_actions[i]);
        });
    }
}

// Methods for submitting the gene form via buttons:
function submitGeneForm(event, form, path) {
    if (titleInput.value === "") {
        if (path === '/form/gene_id/alias') {
            titleInput.value = 'CESA1';
        } else {
            titleInput.value = 'CESA';
        }
    }
    form.submit();
}

function submitNameForm(event, form) {
    const titleInput = document.getElementById("author");
    if (titleInput.value === "") {
        titleInput.value = 'Marek Mutwil';
    }
    form.submit();
}

function submitTitleForm(event, form) {
    const titleInput = document.getElementById("title");
    if (titleInput.value === "") {
        titleInput.value = '26503768';
    }
    form.submit();
}

function submitFormWith(event, value) {
    event.preventDefault();
    // Submit the form with the clicked value
    const link = event.target;
    const form = link.closest('form');
    const input = form.querySelector('input[type="text"]');
    input.value = value;
    form.submit();
}

window.addEventListener('load', () => {
    alias_button_events();
    gene_form_listeners();

    forms[1].setAttribute('onsubmit', `submitNameForm(event, forms[1])`);
    forms[2].setAttribute('onsubmit', `submitTitleForm(event, forms[2])`);
})

/*
This part of the file contains JavaScript for the help part of the landing page.
*/

// Selects the necessary HTML elements to change / attach event listeners to:
const help_information = document.querySelector('#help-text');
const buttons = document.querySelector('.button-group.hollow').querySelectorAll('.button');

// Change help information depending on which button is clicked:
const change_help_text = (text) => {
    let help_text = '';
    switch (text) {
        case 'word':
            help_text = 'Finds all entities containing the searched term.'
            break;
        case 'exact':
            help_text = 'Finds entities containing the exact term (e.g., searching for "cesa" would return "cesa expression", "cesa tracking", and so on).';
            break;
        case 'alias':
            help_text = 'Finds gene aliases (e.g., searching for "cesa1" finds terms such as "rsw1", "atcesa1", and "any1").';
            break;
        case 'substring':
            help_text = 'Finds entities containing the substring (e.g., searching for "hair" would return "root hairs", "hairy roots", and so on).';
            break;
        case 'non-alphanumeric':
            help_text = '"fer" finds entities such as "fer-tor" and "fer/ripk".';
            break;
        default:
            help_text = 'Click one of the above buttons to find out more about each search function!';
            break;
    };
    help_information.innerText = help_text;
}

// Attach event listeners to each button:
buttons.forEach(button => {
    let button_text = button.innerText.toLowerCase();
    button.setAttribute('onclick', `change_help_text('${button_text}')`);
})
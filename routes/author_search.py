'''
This module contains the route(s) needed to search based on an author's name.
'''
from flask import Blueprint, request, render_template
import pickle 
from Bio import Entrez
import sys

# -- Setting up the utils path module -- 
sys.path.append('utils')

# -- Importing custom utilities --
from utils.search import Gene, make_abbreviations, make_functional_annotations
from utils.cytoscape import process_network, generate_cytoscape_js
from utils.text import make_text
from utils.db import fetch_author_ids, fetch_entities_by_ids

author_search = Blueprint('author_search', __name__)

@author_search.route('/author/<query>', methods = ['GET'])
def author(query):
    try:
        my_search = query.lower().strip()
        replacements = {"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss", "é": "e", "ô": "o", "î": "i", "ç": "c"}
        my_search = ''.join(replacements.get(c, c) for c in my_search)
    except:
        my_search = 'Marek Mutwil'.lower()

    if len(my_search):
        hits = fetch_author_ids(' '.join([i.strip() for i in my_search.split()]))
                    
        # provide your email address to the Entrez API
        Entrez.email = "mutwil@gmail.com"

        # search query to find papers by an author
        search_query = my_search + "[Author]"

        # perform the search and retrieve the count of papers
        handle = Entrez.esearch(db = "pubmed", term = search_query)
        record = Entrez.read(handle)
        count = record["Count"]
        
        forSending = []
        if len(hits):
            queries, elements, papers = fetch_entities_by_ids(hits), [], []
            papers = [i[-1] for i in queries]
            for i in queries:
                papers.append(i[-1])
                forSending.append(Gene(i[0], i[2], i[1], i[3])) #source, target, type
                elements.append((i[0].replace("'", "").replace('"', ''),  i[2].replace("'", "").replace('"', ''), i[1].replace("'", "").replace('"', ''))) 

    if len(forSending):
        elements = list(set(elements))
        # fa, ab = pickle.load(open('dbs/fa', 'rb'))[0], pickle.load(open('dbs/abbreviations', 'rb'))[0]
        elementsAb, elementsFa = make_abbreviations(elements), make_functional_annotations(elements)        
        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(updatedElements, elementsAb, elementsFa)
        warning = ''
        summaryText = make_text(forSending)
        
        if len(elements) > 500:
            warning = 'The network might be too large to be displayed, so click on "Layout Options", select the edge types that you are interested in and click "Recalculate layout".'

        return render_template('author.html', genes = forSending, cytoscape_js_code = cytoscape_js_code, 
                               ncbi_count = count, author = query, connectome_count = len(set(papers)), 
                               warning = warning, summary = summaryText)
    else:
        return render_template('not_found.html', search_term = query)

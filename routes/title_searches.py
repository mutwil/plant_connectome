'''
This module contains the route(s) needed to search based on an author's name.
'''
from flask import Blueprint, request, render_template
import pickle
import sys

# -- Setting up the utils path module -- 
sys.path.append('utils')

# -- Importing custom utilities --
from utils.search import Gene, make_abbreviations, make_functional_annotations
from utils.cytoscape import process_network, generate_cytoscape_js
from utils.text import make_text
from utils.db import fetch_title_ids, fetch_entities_by_ids

title_searches = Blueprint('title_searches', __name__)

@title_searches.route('/title/<query>', methods = ['GET'])
def title_search(query):
    try: 
        my_search = query
    except: 
        my_search = '26503768'
    pmids = [i.strip() for i in my_search.split(';')]
            
    forSending = []
    if len(pmids):
        hits, elements = fetch_title_ids(pmids), []
        results = fetch_entities_by_ids(hits)
        for i in results:
            forSending.append(Gene(i[0], i[2], i[1], i[3]))
            elements.append((i[0].replace("'", "").replace('"', ''), i[2].replace("'", "").replace('"', ''), i[1].replace("'", "").replace('"', '')))    

    if forSending!=[]:
        elements = list(set(elements))
        elementsAb, elementsFa = make_abbreviations(elements), make_functional_annotations(elements)

        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(updatedElements, elementsAb, elementsFa)
        summaryText = make_text(forSending)
        return render_template('gene.html', genes = forSending, cytoscape_js_code = cytoscape_js_code, 
                               number_papers = len(hits), search_term = query, summary = summaryText, is_node = False)
    else:
        return render_template('not_found.html', search_term = query)

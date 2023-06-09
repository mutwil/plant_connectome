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

# -- Making Blueprint for title search route:
title_searches = Blueprint('title_searches', __name__)

@title_searches.route('/title/<query>', methods = ['GET'])
def title(query):
    try: 
        my_search = query.lower()
    except: 
        my_search = '26503768'
    pmids = []
    for i in my_search.split(';'):
        pmids += i.split()
            
    forSending = []
    if pmids != []:

        with open('titles', 'rb') as title:
            papers = pickle.load(title)

        hits = list(set(pmids) & set(papers))
        
        
        if hits!=[]:
            with open('allDic2', 'rb') as file:
                genes = pickle.load(file)
            
            
            elements = [] 
            for i in genes:
                for j in genes[i]:
                    if j[3] in hits:
                        if j[0] != '' and j[2] != '':
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type 
                            elements.append((j[0].replace("'", "").replace('"', ''), j[2].replace("'", "").replace('"', ''), j[1].replace("'", "").replace('"', '')))                
                        break
    if forSending!=[]:
        elements = list(set(elements))
        fa, ab = pickle.load(open('fa', 'rb'))[0], pickle.load(open('abbreviations', 'rb'))[0]
        elementsAb, elementsFa = make_abbreviations(ab, elements), make_functional_annotations(fa, elements)

        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(updatedElements, elementsAb, elementsFa)
        summaryText = make_text(forSending)
        return render_template('gene.html', genes = forSending, cytoscape_js_code = cytoscape_js_code, 
                               number_papers = len(hits), search_term = my_search, summary = summaryText)
    else:
        return render_template('not_found.html', search_term = my_search)
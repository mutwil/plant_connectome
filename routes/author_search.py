'''
This module contains the route(s) needed to search based on an author's name.
'''
from flask import render_template, Blueprint
import pickle 
from Bio import Entrez
import sys

# -- Setting up the utils path module -- 
sys.path.append('utils')

# -- Importing custom utilities --
from utils.search import Gene, make_abbreviations, make_functional_annotations
from utils.cytoscape import process_network, generate_cytoscape_js
from utils.text import make_text

author_search = Blueprint('author_search', __name__)

@author_search.route('/author/<query>', methods = ['GET'])
def author(query):
    try:
        my_search = query.lower()
        replacements = {"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss", "é": "e", "ô": "o", "î": "i", "ç": "c"}
        my_search = ''.join(replacements.get(c, c) for c in my_search)
    except:
        my_search = 'Marek Mutwil'.lower()

    if my_search!='':
        with open('authors', 'rb') as f:
            # Load the object from the file
            papers = pickle.load(f)
        hits = []
        for author in papers:      
            if len(set(my_search.split()) & set(author.lower().split())) == len(set(my_search.split())):
                hits += papers[author]
                    
        # provide your email address to the Entrez API
        Entrez.email = "mutwil@gmail.com"

        # search query to find papers by an author
        search_query = my_search + "[Author]"

        # perform the search and retrieve the count of papers
        handle = Entrez.esearch(db = "pubmed", term = search_query)
        record = Entrez.read(handle)
        count = record["Count"]
        
        forSending = []
        if hits != []:
            with open('allDic2', 'rb') as file:
                genes = pickle.load(file)
            
            elements = []
            papers = []
            for i in genes:
                for j in genes[i]:
                    if j[3] in hits:
                        if j[0] != '' and j[2] != '':
                            papers.append(j[3])
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
                            elements.append((j[0].replace("'", "").replace('"', ''),  j[2].replace("'", "").replace('"', ''), j[1].replace("'", "").replace('"', '')))                
                        break
    if forSending != []:
        elements = list(set(elements))
        fa, ab = pickle.load(open('fa', 'rb'))[0], pickle.load(open('abbreviations', 'rb'))[0]
        elementsAb, elementsFa = make_abbreviations(ab, elements), make_functional_annotations(fa, elements)
        
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
        return render_template('not_found.html', search_term = my_search)
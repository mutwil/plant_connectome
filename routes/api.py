from flask import Blueprint, jsonify
from Bio import Entrez
import sys 
import pickle

sys.path.append('utils')
from utils.api import generate_term_api_route, REPLACEMENTS, generate_cytoscape_elements, generate_summary_text
from utils.search import Gene, make_abbreviations, make_functional_annotations
from utils.cytoscape import process_network

api = Blueprint('api', __name__)

# Add these routes in at app.py:
normal_search = generate_term_api_route('default')
exact_search = generate_term_api_route('exact')
substring_search = generate_term_api_route('substring')
alias_search = generate_term_api_route('alias')
nonalpha_search = generate_term_api_route('non-alphanumeric')

# Defining routes for author and title search:
@api.route('/api/author/<query>', methods = ['GET'])
def api_author(query):
	query = ''.join(REPLACEMENTS.get(c, c) for c in query).lower()
	cytoscape_elements, fa, ab = ([], []), [], [],
	text_sum, papers, count, elementsAb, elementsFa = '', 0, 0, [], []

	if len(query):
		papers, hits, forSending = pickle.load(open('dbs/authors', 'rb')), [], []
		for author in papers:
			if len(set(query.split()) & set(author.lower().split())) == len(set(query.split())):
				hits += papers[author]
		
		Entrez.email = 'mutwil@gmail.com'
		search_query = query + "[Author]"
		handle = Entrez.esearch(db = "pubmed", term = search_query)
		record = Entrez.read(handle)
		count = record["Count"]

		if len(hits):
			genes, elements, papers = pickle.load(open('dbs/allDic2', 'rb')), [], []
			for i in genes:
				for j in genes[i]:
					if j[3] in hits and len(j[0]) and len(j[2]):
						papers.append(j[3])
						forSending.append(Gene(j[0], j[2], j[1], j[3]))
						elements.append((j[0].replace("'", "").replace('"', ''),  j[2].replace("'", "").replace('"', ''), j[1].replace("'", "").replace('"', '')))
	
	if len(forSending):
		elements = list(set(elements))
		fa, ab = pickle.load(open('dbs/fa', 'rb'))[0], pickle.load(open('dbs/abbreviations', 'rb'))[0]
		elementsAb, elementsFa = make_abbreviations(ab, elements), make_functional_annotations(fa, elements)       
		cytoscape_elements = generate_cytoscape_elements(process_network(elements))
		text_sum = generate_summary_text(forSending)

	return(jsonify({
		'paper_counts' : {
			'PlantConnectome' : len(set(papers)),
			'NCBI' : int(count)
		},
		'abbreviations' : elementsAb,
		'functional_annotations' : elementsFa,
		'cytoscape_elements' : cytoscape_elements[0] + cytoscape_elements[1],
		'text_summary' : text_sum
	}))

@api.route('/api/title/<query>', methods = ['GET'])
def title_search(query):
	pmids, cytoscape_elements, fa, ab = [i.strip() for i in query.split(';')], ([], []), [], []
	text_sum, elementsAb, elementsFa = '', [], []

	if len(pmids):
		papers = pickle.load(open('dbs/titles', 'rb'))
		hits, forSending = list(set(pmids) & set(papers)), []

		if len(hits):
			genes, elements, papers = pickle.load(open('dbs/allDic2', 'rb')), [], []
			for i in genes:
				for j in genes[i]:
					if j[3] in hits and len(j[0]) and len(j[2]):
						papers.append(j[3])
						forSending.append(Gene(j[0], j[2], j[1], j[3]))
						elements.append((j[0].replace("'", "").replace('"', ''),  j[2].replace("'", "").replace('"', ''), j[1].replace("'", "").replace('"', '')))
	
		if len(forSending):
			elements = list(set(elements))
			fa, ab = pickle.load(open('dbs/fa', 'rb'))[0], pickle.load(open('dbs/abbreviations', 'rb'))[0]
			elementsAb, elementsFa = make_abbreviations(ab, elements), make_functional_annotations(fa, elements)        
			cytoscape_elements = generate_cytoscape_elements(process_network(elements))
			text_sum = generate_summary_text(forSending)

	return(jsonify({
		'abbreviations' : elementsAb,
		'functional_annotations' : elementsFa,
		'cytoscape_elements' : cytoscape_elements[0] + cytoscape_elements[1],
		'text_summary' : text_sum
	}))
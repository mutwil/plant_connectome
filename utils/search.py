'''
This module contains helper utilities to aid in searching (i.e., I've taken YS' and prof.'s code and put 
them in here so that it's easier to work with).
'''
import regex as re
import pickle
from flask import request, render_template

## -- CUSTOM UTILITIES --
from cytoscape import generate_cytoscape_js, process_network
from text import make_text

class Gene:    
    def __init__(self, id, description, inter_type, publication):
        self.id = id
        self.target = description
        self.inter_type = inter_type
        self.publication = publication
    
    def __repr__(self):
        return str(self.__dict__)
    
    def __str__(self):
        return str(self.__dict__)
    
    def getElements(self):
        return (self.id, self.target, self.inter_type)

def is_alphanumeric_helper(string, substring):
    '''
    Checks to see if there are any matches given a string
    and a substring.
    '''
    patternLeft = r'[a-zA-Z0-9]'+ re.escape(substring)
    patternRight = re.escape(substring) + r'[a-zA-Z0-9]'
    matchesLeft = re.findall(patternLeft, string)
    matchesRight = re.findall(patternRight, string)
    
    return (len(matchesLeft) > 0) or (len(matchesRight) > 0)

def find_terms_helper(gene, genes):
    '''
    Converts found terms to a Gene object (defined above)!
    '''
    forSending = []
    elements = []
    for j in genes[gene]:
        if j[0] != '' and j[2] != '':
            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type, pubmed ID
            elements.append((j[0].replace("'", "").replace('"', ''), j[2].replace("'", "").replace('"', ''), 
                             j[1].replace("'", "").replace('"', '')))
    return elements, forSending

def find_terms(my_search, genes, search_type):   
    '''
    Given a search term, something to search (i.e., a pickled Python file), 
    and a search type, return the items that match the search query.

    This function is KEY to the searching functionality on the application!
    '''
    if not len(my_search):
        return None
    
    forSending, elements = [], []
    if search_type == 'exact':
        for i in genes:
            if my_search.upper().strip() == i.strip():
                elements, forSending = find_terms_helper(i, genes)
    elif search_type == 'alias':
        with open('dbs/geneAlias', 'rb') as file:
            adjMat = pickle.load(file)
        try:
            terms = adjMat[my_search.upper().strip()]
        except:
            terms = []
        for i in terms:
            for j in genes:
                if i.upper().strip() in j.strip().split():
                    outputOne, outputTwo = find_terms_helper(j, genes)
                    elements.extend(outputOne)
                    forSending.extend(outputTwo)
    elif search_type == 'substring':
        for i in genes:
            if my_search.upper().strip() in i.strip():
                outputOne, outputTwo = find_terms_helper(i, genes)
                elements.extend(outputOne)
                forSending.extend(outputTwo)
    elif search_type == 'non-alphanumeric':
        for i in genes:
            substring = my_search.upper().strip()
            string = i.strip()
            if substring in string:
                if not is_alphanumeric_helper(string, substring):
                    outputOne, outputTwo = find_terms_helper(i, genes)
                    elements.extend(outputOne)
                    forSending.extend(outputTwo)
    elif search_type == 'default':                                                       # (i.e., default case)
        for i in genes:
            if my_search.upper().strip() in i.strip().split():  # default search - terms that contain the specific query word
                outputOne, outputTwo = find_terms_helper(i, genes)
                elements.extend(outputOne)
                forSending.extend(outputTwo)
    else:
        raise Exception("An invalid 'search_type' parameter has been entered!")
    return list(set(elements)), forSending

def make_abbreviations(abbreviations, elements):
    ab = {}
    for element in elements:
        if abbreviations.get(element[0].upper()) is not None:
            ab[element[0]] = abbreviations[element[0].upper()]
        if abbreviations.get(element[1].upper()) is not None:
            ab[element[1]] = abbreviations[element[1].upper()] 
    return ab

def make_functional_annotations(gopredict, elements):
    fa = {}
    for element in elements:
        if gopredict.get(element[0].upper()) is not None:
            fa[element[0]] = gopredict[element[0].upper()]
        if gopredict.get(element[1].upper()) is not None:
            fa[element[1]] = gopredict[element[1].upper()]
    return fa
    
def generate_search_route(search_type):
    '''
    A function factory that generates a route to be used in the flask application.  
    The routes for the search form for genes, aliases, and other entities are all 
    virtually similar to one another - a function factory keeps things DRY.
    '''
    def search_route(query):
        try:
            my_search = query
        except:
            my_search = 'cesa'
        
        forSending = []
        if len(my_search) > 0:
            split_search = my_search.split(';')
            elements = []
  
            to_search = pickle.load(open('dbs/allDic2', 'rb'))
            ab = pickle.load(open('dbs/abbreviations', 'rb'))[0]
            fa = pickle.load(open('dbs/fa', 'rb'))[0]

            for term in split_search:
                results = find_terms(term, to_search, search_type)
                elements += results[0]
                forSending += results[1]
                
            # remove redundancies
            elements = list(set(elements))
            
            warning = ''
            if len(elements) > 500:
                warning = 'The network might be too large to be displayed, so click on "Layout Options",  select the edge types that you are interested in and click "Recalculate layout".'
            
            updatedElements = process_network(elements)
            elementsAb = make_abbreviations(ab, elements)
            elementsFa = make_functional_annotations(fa, elements)
            cytoscape_js_code = generate_cytoscape_js(updatedElements, elementsAb, elementsFa)
            if elementsAb.get(my_search.upper()) is not None:
                node_ab = elementsAb[my_search.upper()]
            else:
                node_ab = []

            if elementsFa.get(my_search.upper()) is not None:
                node_fa = elementsFa[my_search.upper()]
            else:
                node_fa = []
            
            papers = []
            for i in forSending:
                papers += [i.publication]
                
            summaryText = make_text(forSending)
            
        if forSending != []:
            return render_template('gene.html', genes = forSending, cytoscape_js_code = cytoscape_js_code, 
                                    search_term = query, number_papers = len(set(papers)), warning = warning, 
                                    summary = summaryText, node_ab = node_ab, node_fa = node_fa, is_node = True)
        else:
            return render_template('not_found.html', search_term = query)
    return search_route

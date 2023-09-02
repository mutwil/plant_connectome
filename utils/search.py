'''
This module contains helper utilities to aid in searching (i.e., I've taken YS' and prof.'s code and put 
them in here so that it's easier to work with).
'''
import regex as re
from flask import request, render_template
from pymongo import MongoClient

## -- CUSTOM UTILITIES --
from cytoscape import generate_cytoscape_js, process_network
from text import make_text

URL = 'mongodb://localhost:27017/'
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

def find_terms_helper(searches):
    '''
    Converts found terms to a Gene object (defined above)!
    '''
    forSending, elements = [], []
    for i in searches:
        if len(i['entity1']) and len(i['entity2']):
            forSending.append(Gene(i['entity1'], i['entity2'], i['edge'], i['pubmedID']))
            elements.append((i['entity1'].replace("'", ''), i['entity2'].replace("'", ''), 
                             i['edge'].replace("'", ''), i['pubmedID'].replace("'", '')))
    return elements, forSending

def find_terms(my_search, search_type, url = 'mongodb://localhost:27017/'):   
    '''
    Given a search term, something to search (i.e., a pickled Python file), 
    and a search type, return the items that match the search query.

    This function is KEY to the searching functionality on the application!
    '''
    client = MongoClient(url)
    table = client.plant_connectome['all_dic']
    if not len(my_search):
        return None
    
    forSending, elements = [], []
    if search_type == 'exact':
        results = table.find({'$or' : [{'entity1' : my_search.upper()}, {'entity2' : my_search.upper()}]})
        elements, forSending = find_terms_helper(results)
    elif search_type == 'alias':
        alias_results, terms = list(client.plant_connectome['gene_alias'].find({'gene' : my_search.upper()})), []
        for i in alias_results:
            terms.extend(i['aliases']) 
        for i in terms:
            results = table.find({"$or" : [{'entity1' : {'$regex' : fr'\b{i.upper()}\b'}}, 
                                       {'entity2' : {'$regex' : fr'\b{i.upper()}\b'}}]})
            outputOne, outputTwo = find_terms_helper(results)
            elements.extend(outputOne) ; forSending.extend(outputTwo)
    elif search_type == 'substring':
        results = table.find({"$or" : [{'entity1' : {'$regex' : fr'\b{my_search.upper()}\b'}}, 
                                       {'entity2' : {'$regex' : fr'\b{my_search.upper()}\b'}}]})
        elements, forSending = find_terms_helper(results)
    elif search_type == 'non-alphanumeric':
        results = list(table.find({"$or" : [{'entity1' : {'$regex' : fr'\b{my_search.upper()}\b'}}, 
                                       {'entity2' : {'$regex' : fr'\b{my_search.upper()}\b'}}]}))
        results_copy = results
        for i in results_copy: 
            if not is_alphanumeric_helper(i['entity1'], my_search.upper().strip()):
                continue
            if not is_alphanumeric_helper(i['entity2'], my_search.upper().strip()):
                continue
            results_copy.append(i)
        elements, forSending = find_terms_helper(results_copy)
    elif search_type == 'default':                                                       # (i.e., default case)
        results = table.find({"$or" : [{'entity1' : {'$regex' : fr'\b{my_search.upper()}\b'}}, 
                                       {'entity2' : {'$regex' : fr'\b{my_search.upper()}\b'}}]})
        elements, forSending = find_terms_helper(results)
    else:
        raise Exception("An invalid 'search_type' parameter has been entered!")
    
    client.close()
    return list(set(elements)), forSending

def make_abbreviations(elements, url = 'mongodb://localhost:27017/'):
    client = MongoClient(url)
    table = client.plant_connectome.abbr
    results = [table.find({'$or' : [{'term' : i[0].upper()}, {'term' : i[1].upper()}]}) for 
              i in elements]
    ab = {j['term'] : j['abbs'] for i in results for j in i}
    client.close()
    return ab

def make_functional_annotations(elements, url = 'mongodb://localhost:27017/'):
    client = MongoClient(url)
    table = client.plant_connectome.fa
    results = [table.find({'$or' : [{'term' : i[0].upper()}, {'term' : i[1].upper()}]}) for 
              i in elements]
    fa = {j['term'] : j['abbs'] for i in results for j in i}
    client.close()
    return fa
    
def generate_search_route(search_type):
    '''
    A function factory that generates a route to be used in the flask application.  
    The routes for the search form for genes, aliases, and other entities are all 
    virtually similar to one another - a function factory keeps things DRY.
    '''
    def search_route(query):
        forSending, my_search = [], query

        if len(my_search) > 0:
            split_search, elements = my_search.split(';'), []
            for term in split_search:
                results = find_terms(term, search_type)
                elements += results[0] ; forSending += results[1]
                
            elements = list(set(elements))
            
            warning = ''
            if len(elements) > 500:
                warning = 'The network might be too large to be fully displayed.'
            
            updatedElements = process_network(elements)
            elementsAb, elementsFa = make_abbreviations(elements), make_functional_annotations(elements)
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

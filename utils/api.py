'''
Contains API-related information that the api route use...
'''
from flask import jsonify
from search import find_terms, make_abbreviations, make_functional_annotations
from cytoscape import process_network
import pickle

REPLACEMENTS = {"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss", "é": "e", "ô": "o", "î": "i", "ç": "c"}

def generate_cytoscape_elements(elements):
    '''
    Generates two lists that can be JSON-ed 
    and put into a CytoscapeJS session.
    '''
    nodes, edges = [], []
    for node in set(edge["source"] for edge in elements) | set(edge["target"] for edge in elements):
        nodes.append({'data' : {'id' : node}})

    for i, edge in enumerate(elements):
        edges.append({
            'data' : {'id' : i, 'source' : edge['source'], 'target' : edge['target'], 'interaction' : edge['interaction']}
        })
    return nodes, edges

def generate_summary_text(elements):
    if not len(elements):
        return ''
    
    topicDic, nodeDegree, nodeSentenceDegree = {}, {}, {}
    for i in elements:
        if i.id not in topicDic:
            topicDic[i.id] = {}
            nodeDegree[i.id] = 0
            nodeSentenceDegree[i.id] = {}
        if i.inter_type not in topicDic[i.id]:
            topicDic[i.id][i.inter_type] = [[i.target, i.publication]]
            nodeSentenceDegree[i.id][i.inter_type] = 1
            nodeDegree[i.id] += 1
        else:
            topicDic[i.id][i.inter_type] += [[i.target, i.publication]]
            nodeDegree[i.id] += 1
            nodeSentenceDegree[i.id][i.inter_type] += 1

    sorted_nodes, finished_sentences = sorted(nodeDegree, key = lambda x: nodeDegree[x], reverse = True), []
    for i in sorted_nodes:
        sorted_sentences = sorted(nodeSentenceDegree[i], key = lambda x: nodeSentenceDegree[i][x], reverse = True)
        for j in sorted_sentences:
            text, temp, tempRefs = i + ' ' + j + ' ', {}, []
            for k in topicDic[i][j]:
                if k[0] not in temp:
                    temp[k[0]] = [k[1]]
                else:
                    temp[k[0]] += [k[1]]
            
            for target in temp:
                tempRefs += [target + ' (' + ', '.join(list(set(temp[target])))+ ')']
            finished_sentences.append(text + ', '.join(tempRefs))

    return '. '.join(finished_sentences) + '.'


def generate_term_api_route(query_type):
    '''
    Returns a function meant to be used a view handler in the APIs.
    '''
    def api_route(query):
        if len(query):
            forSending, elements, summary = [], [], ''
            to_search, func_anno = pickle.load(open('dbs/allDic2', 'rb')), pickle.load(open('dbs/fa', 'rb'))[0]
            abbr = pickle.load(open('dbs/abbreviations', 'rb'))[0]

            for term in query.split(';'):
                results = find_terms(term, to_search, query_type)
                elements.extend(results[0]) ; forSending.extend(results[1])

            elements, papers, summary = list(set(elements)), [i.publication for i in forSending], generate_summary_text(forSending)
            elementsAb, elementsFa = make_abbreviations(abbr, elements), make_functional_annotations(func_anno, elements)
            cytoscape_elements = generate_cytoscape_elements(process_network(elements))
        
        return jsonify({
            'functional_annotations' : elementsFa,
            'abbreviations' : elementsAb,
            'cytoscape_elements' : cytoscape_elements[0] + cytoscape_elements[1],
            'text_summary' : summary,
            'publications' : papers
        })
        
    return api_route
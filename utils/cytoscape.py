'''
This module contains helper functions needed to generate the CytoscapeJS graph.
'''
import networkx as nx
import json

def graphConverter(graph, ref):
    updatedElements = []
    for k, v in graph.adjacency():
        source = str(k[0]).replace("'", "").replace('"', '')
        if v:
            for i, j in v.items():
                target = str(i).replace("'", "").replace('"', '')
                for p, q in j.items():
                    type = str(q['relation']).replace("'", "").replace('"', '')
                    if ((source, k[1]), target) in ref:
                        updatedElements.append({"source": source, "target": target, "interaction": type})
    return updatedElements 
            
def edgeConverter(elements): #Convert Edges to default dictionary format 
    updatedElements = []
    for i in elements:
        updatedElements.append({"source": str(i[0]).replace("'", "").replace('"', ''), 
                                "target": str(i[1]).replace("'", "").replace('"', ''), 
                                "interaction": str(i[2]).replace("'", "").replace('"', '')})
    return updatedElements 

def nodeDegreeFilter(graph):
    nodesToKeep = []
    totalDegree = 0
    for i in sorted(dict(graph.degree), key = dict(graph.degree).get, reverse = True):
        currDegree = graph.degree(i)
        if totalDegree <= 500 or nodesToKeep == []: #must keep at least one nodes
            nodesToKeep.append(i)
            totalDegree += currDegree
        if totalDegree > 500:
            break
    lst = set()
    for j in nodesToKeep:
        lst.update(list(graph.in_edges(j)) + list(graph.out_edges(j)))
    ref = {}
    for i in lst:  # create a hash map for the entities to be kept
        ref[i] = ref.get(i, 0) + 1
    return graph, ref

def process_network(elements):
    elements = list(set(elements))    
    if len(elements) > 500:  # only perform node filtration if the number of relationship is larger than 500
        #using Networkx multiDiGraph to model input data as Graph
        G = nx.MultiDiGraph()
        for i in elements:
            G.add_edge((i[0], i[-1]), i[1], relation = i[2])
        
        G, ref = nodeDegreeFilter(G)
        
        return graphConverter(G, ref)
    else:
        return edgeConverter(elements)

def generate_cytoscape_js(elements, ab, fa):
    '''
    Generates nodes and edges to be displayed in the Cytoscape JS network.
    '''
    nodes = [
        "{ data: { id: '%s' } }" % node
        for node in set(edge["source"] for edge in elements) | set(edge["target"] for edge in elements)
    ]

    edges = [
        "{ data: { id: 'edge%s', source: '%s', target: '%s', interaction: '%s' } }" % (
            i,
            edge['source'],
            edge['target'],
            edge['interaction'],
        )
        for i, edge in enumerate(elements)
    ]
    
    a = open('stats/network.txt', 'r').read()

    nodes = ', '.join(nodes)
    edges = ', '.join(edges)
    
    return a.replace('NODES',nodes).replace('EDGES',edges).replace('REPLACE_AB', json.dumps(ab)).replace('REPLACE_FA', json.dumps(fa))

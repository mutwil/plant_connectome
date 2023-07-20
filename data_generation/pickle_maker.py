# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 21:58:11 2023

@author: marek
"""
import os
import collections
import pickle

def simplify_edge(edge: str, to_search: dict) -> str:
    '''
    Determines the simplified edge that an edge should be transformed into given its value.  This function is 
    heavily dependent on the SIMPLIFICATIONS dictionary inside edges.py, so as time goes on, do add more groupings and their associated 
    edges!

    If no grouping is found, then the original edge is returned
    '''
    for pairs in list(to_search.items()):
        if edge.upper() in pairs[1]:
            return pairs[0].upper()
    return edge.upper()

def simplify_passive(edge: str, to_search: dict) -> str:
    '''
    Determines the simplified edge that an edge should be transformed into given its value.  This function is 
    heavily dependent on the SIMPLIFICATIONS dictionary inside edges.py, so as time goes on, do add more groupings and their associated 
    edges!

    If no grouping is found, then the original edge is returned
    '''
    for pairs in list(to_search.items()):
        if edge.upper() in pairs[1]:
            return pairs[0].upper()
    return ''

# path = 'D:/GPT_annotator/'

# GROUPINGS = edges.GROUPINGS
merges = pickle.load(open('merged_edges', 'rb'))
passive_merges = pickle.load(open('passive_edges', 'rb'))

'''
BASICALLY... it appends the edges (i.e., the text in between two exclamation marks) to 
the list 'types'.
'''
types = []
for i in os.listdir('../annotations/'):
    try:
        a = open('../annotations/' + i).read()
        pmid, abstract, pubdate, journal, doi, main = eval(a.split('\n')[0])
        findings = a.split('\n\n')[1].split(' \n')
        for j in findings:
            if j.count('!') == 2:
                types.append(j.split('!')[1])   # Appends edges
    except:
        pass
'''
interactions_count = collections.Counter(types)
top_studied = dict(interactions_count.most_common(100))

I made interaction_types from top 100 types
'''

'''
types2type = {}

for i in open('D:/GPT_annotator/interaction_types.txt','r').readlines():
    if ':' not in i:
        types2type[i.rstrip()] = i.rstrip()
    else:
        splitta = i.rstrip().split(':')
        types2type[splitta[0].rstrip()] = splitta[0].rstrip()
        for j in splitta[1].split(','):
            types2type[j.strip()] = splitta[0].rstrip()
'''

# === Makign the allDic pickle file ===
'''
UPDATE 05/05/2023 --> prof. Marek suggests grouping the terms here:
'''

good, bad = 0,0

agent2agents = {}
uniqueEdges = []
errors = set()
for i in os.listdir('../annotations/'):
    try:
        a = open('../annotations/' + i, encoding = 'iso-8859-1').read()
        
        findings = a.split('\n\n')[1].split('\n')
        
        for j in findings:
            if j.count('!') == 2:
                splitta = j.split('!')
    
                agentA, _, agentB = splitta
                agentA = agentA.split(':')[0].strip().upper()
                agentB = agentB.strip().upper()
                edge, active_edge = simplify_edge(_.strip(), merges), simplify_passive(_.strip(), passive_merges)
                '''
                A tuple that contains the following information:
                (node1, edge, node2, name of the file)
                '''
                if len(agentB) > 1 and len(agentA) > 1:
                    if len(active_edge):
                        edge_tuple = (agentB.upper(), active_edge, agentA.upper(), i.split('.')[0])
                    else: 
                        edge_tuple = (agentA.upper(), edge, agentB.upper(), i.split('.')[0])
                else:
                    print(f"SKIPPING THE FOLLOWING PAIR (too short): {agentA if not len(agentA) else '<EMPTY>'} and {agentB if not len(agentB) else '<EMPTY>'}!")
                    continue 

                #if agent absent
                if agentA not in agent2agents:
                    agent2agents[agentA] = []
                if agentB not in agent2agents:
                    agent2agents[agentB] = []
                
                #the problem is here: rsw1 synthesized cellulose, cellulose synthesized rsw1
                agent2agents[agentA] += [edge_tuple]
                agent2agents[agentB] += [edge_tuple]
                uniqueEdges.append(edge_tuple)
    except Exception as e:
        errors.add(str(e))
        print(i, 'wonky')

print(errors)

with open('../dbs/allDic2', 'wb') as file:
    pickle.dump(agent2agents, file)

# =====

save = []
for i in list(set(uniqueEdges)):
    save.append(i[0]+'\t'+i[2]+'\n')
  
try:
    with open('edges.txt', 'w', encoding='utf-8') as f:
        for string in save:
            f.write(string + '\n')
except:
    print('count not save')
    

items, edges = [],0
for i in os.listdir('../annotations/'):
    a = open('../annotations/'+i,'r').read()
    
    if len(a)>0:
        findings = a.split('\n\n')[1].split('\n')
        
        for j in findings:
            if j.count('!')==2:
                splitta = j.split('!')
    
                agentA, _, agentB = splitta
                agentA = agentA.split(':')[0].upper()
                agentB = agentB.strip().upper()
                edges+=1
                items+=[agentA]
                items+=[agentB]
            

v = open('../stats/stats.txt','w')
v.write(str(len(os.listdir('../annotations/')))+'\t'+str(len(set(items))))
v.close()


papers = []
with open('../dbs/abstracts', 'rb') as f:
    # Load the object from the file
    abstracts = pickle.load(f)
    for j in abstracts:
        papers.append(j['journal'])
        
counted_items = collections.Counter(papers)
top_20_items = counted_items.most_common(20)

keys = []
values = []

for item in top_20_items:
    keys.append(item[0])
    values.append(item[1])

keys.append('Other')
values.append(len(papers)-sum(values))
with open('../stats/journal_statistics.txt', 'w') as f:
    f.write(str(keys) + '\n')
    f.write(str(values) + '\n')

save = {}
replacements = {"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss", "é": "e", "ô": "o", "î": "i", "ç": "c"}
for paper in abstracts:
    for j in paper['authors']:            
        author = ''.join(replacements.get(c, c) for c in j)
        if j not in save:
            save[j] = []
        save[j]+=[paper['pmid']]
        
with open('../dbs/authors', 'wb') as file:
    pickle.dump(save, file)
    
    
save = []
for paper in abstracts:
    save.append(paper['pmid'])
        
with open('../dbs/titles', 'wb') as file:
    pickle.dump(set(save), file)
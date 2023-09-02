'''
A module that includes utility functions to 
assist with fetching information from the MongoDB database.

== NOT IN USE NOW ==
'''
from pymongo import MongoClient

def fetch_entities(con = 'mongodb://localhost:27017', db = 'plant_connectome'):
    '''
    Returns all entities in the cata table.
    '''
    client = MongoClient(con)
    cata = list(client[db].cata.find({}))
    letters, items = [i['letter'] for i in cata], {i['letter'] : i['terms'] for i in cata}
    client.close()
    return letters, items

def fetch_author_ids(author, con = 'mongodb://localhost:27017', db = 'plant_connectome'):
    '''
    Returns all ids in the authors table that match the
    'author' search_query.
    '''
    client = MongoClient(con)
    result, ids = list(client[db].authors.find({'author' : author.lower()})), []
    client.close()
    for i in result:
        ids.extend(i['ids'])
    return ids

def fetch_entities_by_ids(ids, con = 'mongodb://localhost:27017', db = 'plant_connectome'):
    client = MongoClient(con)
    query = list(client[db].all_dic.find({'pubmedID' : {'$in' : ids}}))
    client.close()
    return [(i['entity1'], i['entity2'], i['edge'], i['pubmedID']) for i in query if len(i['entity1']) and len(i['entity2'])]

def fetch_title_ids(ids, con = 'mongodb://localhost:27017', db = 'plant_connectome'):
    client = MongoClient(con)
    results = client[db].titles.find({'id' : {'$in' : list(set(ids))}})
    client.close()
    return [i['id'] for i in results]

def find_abs(query, con = 'mongodb://localhost:27017', db = 'plant_connectome'):
    client = MongoClient(con)
    query = client[db].abbr.find({'term' : query.upper()})
    client.close()
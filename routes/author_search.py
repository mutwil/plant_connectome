'''
This module contains the route(s) needed to search based on an author's name.
'''
from flask import Blueprint, request, render_template
import pickle
from Bio import Entrez
import sys

# -- Setting up the utils path module --
sys.path.append('utils')
from utils.mongo import db
from utils.text import make_text
from utils.search import Gene, make_abbreviations, make_functional_annotations
from utils.cytoscape import process_network, generate_cytoscape_js
# -- Importing custom utilities --


author_search = Blueprint('author_search', __name__)


@author_search.route('/author/<query>', methods=['GET'])
def author(query):
    try:
        my_search = query.lower()
        replacements = {"ä": "ae", "ö": "oe", "ü": "ue",
                        "ß": "ss", "é": "e", "ô": "o", "î": "i", "ç": "c"}
        my_search = ''.join(replacements.get(c, c) for c in my_search)
    except:
        my_search = 'Mutwil Marek'.lower()
# mongo integration =>
    authors_collection = db["authors"]
    all_dic_collection = db["all_dic"]
    # fa_collection = db["fa"]
    # ab_collection = db["abbr"]
    # all_view = db.all_view
    forSending = []
    elements = []
    elementsAb = {}
    elementsFa = {}
    if my_search != '':
        hits = authors_collection.find_one({"author": my_search})
        # provide your email address to the Entrez API
        Entrez.email = "mutwil@gmail.com"

        # search query to find papers by an author
        search_query = my_search + "[Author]"

        # perform the search and retrieve the count of papers
        handle = Entrez.esearch(db="pubmed", term=search_query)
        record = Entrez.read(handle)
        count = record["Count"]
        if hits is not None:
            result = all_dic_collection.find(
                {"pubmedID": {"$in": hits["ids"]}})
            for i in result:
                forSending.append(
                    Gene(i["entity1"], i["entity2"], i["edge"], i["pubmedID"]))
                elements.append((i["entity1"],
                                i["entity2"], i["edge"]))
                if elementsAb.get(i["entity1"]) is None and len(i["entity1_abbs"]) > 0:
                    elementsAb[i["entity1"]] = i["entity1_abbs"]
                if elementsFa.get(i["entity1"]) is None and len(i["entity1_fas"]) > 0:
                    elementsFa[i["entity1"]] = i["entity1_fas"]

                if elementsAb.get(i["entity2"]) is None and len(i["entity2_abbs"]) > 0:
                    elementsAb[i["entity2"]] = i["entity2_abbs"]
                if elementsFa.get(i["entity2"]) is None and len(i["entity2_fas"]) > 0:
                    elementsFa[i["entity2"]] = i["entity2_fas"]
# <= mongo integration

        #     with open('authors', 'rb') as f:
        #         # Load the object from the file
        #         papers = pickle.load(f)
        #     hits = []
        #     for author in papers:
        #         if len(set(my_search.split()) & set(author.lower().split())) == len(set(my_search.split())):
        #             hits += papers[author]

            # # provide your email address to the Entrez API
            # Entrez.email = "mutwil@gmail.com"

            # # search query to find papers by an author
            # search_query = my_search + "[Author]"

            # # perform the search and retrieve the count of papers
            # handle = Entrez.esearch(db = "pubmed", term = search_query)
            # record = Entrez.read(handle)
            # count = record["Count"]

        #     forSending = []
        #     if hits != []:
        #         with open('allDic2', 'rb') as file:
        #             genes = pickle.load(file)

        #         elements = []
        #         papers = []
        #         for i in genes:
        #             for j in genes[i]:
        #                 if j[3] in hits:
        #                     if j[0] != '' and j[2] != '':
        #                         papers.append(j[3])
        #                         forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
        #                         elements.append((j[0].replace("'", "").replace('"', ''),  j[2].replace("'", "").replace('"', ''), j[1].replace("'", "").replace('"', '')))
        #                     break
        # if forSending != []:
        #     elements = list(set(elements))
        #     fa, ab = pickle.load(open('fa', 'rb'))[0], pickle.load(open('abbreviations', 'rb'))[0]
        #     elementsAb, elementsFa = make_abbreviations(ab, elements), make_functional_annotations(fa, elements)
        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(
            updatedElements, elementsAb, elementsFa)
        warning = ''
        summaryText = make_text(forSending)

        if len(elements) > 500:
            warning = 'The network might be too large to be displayed, so click on "Layout Options", select the edge types that you are interested in and click "Recalculate layout".'

        return render_template('author.html', genes=forSending, cytoscape_js_code=cytoscape_js_code,
                               ncbi_count=count, author=query, connectome_count=len(
                                   set(hits["ids"])),
                               warning=warning, summary=summaryText)
    else:
        return render_template('not_found.html', search_term=query)

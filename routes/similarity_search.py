'''
This module contains the routes for similarity search.
'''
from flask import Blueprint, request, render_template
import pickle
import sys 

## -- CUSTOM UTILITIES --

# -- Setting up the utils path module -- 
sys.path.append('utils')

similarity_search = Blueprint('similarity_search', __name__)

@similarity_search.route('/similarity', methods = ['POST'])
def similarity():
    try:
        my_search = request.form["similarity_id"]
        type = request.form["similarity_type"]
    except:
        my_search = ""
        type = ""

    forSending = []
    if my_search != "" and type != "":
        # search all the nodes with the same term
        genes = pickle.load(open('allDic2', 'rb'))
        if type == "ab":
            f = pickle.load(open('abbreviations', 'rb'))
            typa = "abbreviation"
        else:
            f = pickle.load(open('fa', 'rb'))
            typa = "functional annotation"

        for k, v in f[0].items():
            if my_search in v and k in f[1]:
                forSending.append((k, f[1][k]))

        unique_papers = []
        for i in forSending:
            for v in i[1]:
                if v not in unique_papers:
                    unique_papers.append(v)
       
        return render_template('similarity.html', results = forSending, search_term = my_search, number_nodes = len(forSending), number_papers = len(unique_papers), type = typa)


    else:
        return render_template('not_found.html', search_term = my_search)


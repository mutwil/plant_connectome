'''
This module contains the routes for entities.
'''
from flask import Blueprint, request, render_template
import pickle
import sys 
import re

## -- CUSTOM UTILITIES --

# -- Setting up the utils path module -- 
sys.path.append('utils')

catalogue_search = Blueprint('catalogue_search', __name__)
@catalogue_search.route('/catalogue', methods = ['GET'])
def catalogue():
    cata = pickle.load(open('dbs/cata2','rb'))
    return render_template("/catalogue.html", entities = cata[1], header=sorted(cata[0]))

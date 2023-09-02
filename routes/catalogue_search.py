'''
This module contains the routes for entities.
'''
from flask import Blueprint, render_template

# -- Setting up the utils path module -- 
from utils.db import fetch_entities

catalogue_search = Blueprint('catalogue_search', __name__)
@catalogue_search.route('/catalogue', methods = ['GET'])
def catalogue():
    letters, entities = fetch_entities()
    return render_template("/catalogue.html", entities = entities, header = sorted(letters))

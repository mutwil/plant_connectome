'''
This module contains the routes for searching.
'''
# from flask import Blueprint
# import pickle
import sys 

# -- Setting up the utils path module -- 
sys.path.append('utils')

# -- Importing custom utilities --
from utils.search import generate_search_route

# term_searches = Blueprint('term_searches', __name__)

# -- Constants --
#DATABASE = pickle.load(open('allDic2', 'rb'))
#ABBREVIATIONS = pickle.load(open('abbreviations', 'rb'))
#FUNCANNOTATE = pickle.load(open('fa', 'rb'))

# -- Generating the routes via function factories:
normal = generate_search_route('default')
exact = generate_search_route('exact')
alias = generate_search_route('alias')
substring = generate_search_route('substring')
non_alpha = generate_search_route('non-alphanumeric')

# -- Adding the routes (commented out for now) --
'''
term_searches.add_url_rule('/normal/<query>', endpoint = 'normal', view_func = normal, methods = ['POST'])
term_searches.add_url_rule('/exact/<query>', endpoint = 'exact', view_func = exact, methods = ['POST'])
term_searches.add_url_rule('/alias/<query>', endpoint = 'alias', view_func = alias, methods = ['POST'])
term_searches.add_url_rule('/substring/<query>', endpoint = 'substring', view_func = substring, methods = ['POST'])
term_searches.add_url_rule('/non_alpha/<query>', endpoint = 'non_alpha', view_func = non_alpha, methods = ['POST'])
'''
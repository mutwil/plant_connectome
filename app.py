from flask import Flask, render_template, request, url_for, redirect, send_from_directory
import os

# == IMPORTING THE ROUTES  ==
from routes.author_search import author
from routes.term_searches import normal, exact, alias, substring, non_alpha
from routes.title_searches import title_searches
from routes.similarity_search import similarity_search
from routes.author_search import author_search

app = Flask(__name__)

# == REGISTERING THE ROUTES ==
app.register_blueprint(author_search)
app.register_blueprint(title_searches)
app.register_blueprint(similarity_search)

# -- REGISTERING FORM SEARCHES AS DYNAMIC ROUTES --
with app.app_context():
    app.add_url_rule('/normal/<query>', 'normal', normal, methods = ['GET'])
    app.add_url_rule('/exact/<query>', 'exact', exact, methods = ['GET'])
    app.add_url_rule('/alias/<query>', 'alias', alias, methods = ['GET'])
    app.add_url_rule('/substring/<query>', 'substring', substring, methods = ['GET'])
    app.add_url_rule('/non_alpha/<query>', 'non_alpha', non_alpha, methods = ['GET'])

@app.route('/form/<form_type>/<search_type>', methods = ['POST'])
def form(form_type, search_type):
    query = request.form[form_type]
    return redirect(url_for(search_type, query = query))

# == Fetching the 'index', the 'help', and the 'features' page ==
'''
I left these routes here as they're trivial enough for now.
'''

@app.route('/', methods = ['GET'])
def index():
    v = open('stats.txt','r').read().rstrip().split()
    return render_template('index.html', entities = v[1], papers = v[0])

@app.route('/help', methods = ['GET'])
def help():
    '''
    Returns the help page.
    '''
    return render_template('help.html')

@app.route('/features', methods = ['GET'])
def features():
    '''
    Renders the features page.
    '''
    journals, numbers = open('journal_statistics.txt','r').read().splitlines()
    piechart = open('piechart.txt','r').read().replace('JOURNALS', journals).replace('NUMBERS', numbers)
    return render_template('features.html', piechart_code = piechart)

# Fetches the favicon:
@app.route('/favicon.ico')
def favicion():
    return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico', mimetype = 'image/vnd.microsoft.icon')

if __name__ == '__main__':
    app.run()
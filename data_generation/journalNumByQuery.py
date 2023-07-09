#@markdown Run this cell to load in get_papers(search_term) function
from Bio import Entrez

import os
import openai
import json
import time
import pickle 

prompt = '''
Write a very short summary about the functions of genes in this abstract. The summary must show pair-wise relationships, for example:
gene: !affects! Process
gene: !localizes to! X
gene: !interacts with! Y
gene: !enhances! Z
gene: !represses! U
gene: !synthesizes! I

Please provide only one statement per line, and ensure that each line contains exactly two actors. If a relationship involves more than two actors, please break it down into multiple separate lines.

<ENTER PROMPT HERE>

VERY SHORT, CONCISE SUMMARY CONTAINING ALL INFORMATION WITH TWO ACTORS PER LINE: 

'''

openai.api_key = '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'
path = 'D:/GPT_annotator'

def gpt_analyzer(abstract):
    print('Abstract:',  abstract)
    my_prompt = prompt.replace('<ENTER PROMPT HERE>',abstract)

    response = openai.Completion.create(
      model="text-davinci-003",
      prompt=my_prompt,
      temperature=0,
      max_tokens=256,
      top_p=1.0,
      frequency_penalty=0.0,
      presence_penalty=1
    )
    print(response['choices'][0]['text']+'\n')

    return response['choices'][0]['text']

def search_and_fetch_journal_papers(query_term, email, year_cutoff, retmax=100000):
    Entrez.email = email
    query = f'"{query_term}" AND "{year_cutoff}/01/01"[Date - Publication] : "3000"[Date - Publication]'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=retmax,
                            term=query,
                            usehistory='y')
    results = Entrez.read(handle)
    handle.close()

    pmid_list = results['IdList']
    handle = Entrez.efetch(db='pubmed', id=','.join(pmid_list), retmode='xml')
    papers = Entrez.read(handle)['PubmedArticle']
    handle.close()

    paper_details = []
    print('pmids', len(papers), query_term)
    for paper in papers:
        pmid = paper['MedlineCitation']['PMID']
        title = paper['MedlineCitation']['Article']['ArticleTitle']
        abstract = paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', None)
        abstract = ' '.join(abstract) if abstract is not None else ''
        date = paper['MedlineCitation']['Article']['ArticleDate'][0] if paper['MedlineCitation']['Article'].get('ArticleDate') else paper['MedlineCitation']['DateRevised']
        journal = query_term.split('[')[0]
        doi = paper['PubmedData']['ArticleIdList'][-1] if paper['PubmedData'].get('ArticleIdList') else None
        author_list = paper['MedlineCitation']['Article'].get('AuthorList', [])
        authors = [f"{author.get('LastName', '')} {author.get('ForeName', '')}" for author in author_list if author.get('LastName')]

        paper_info = {
            'pmid': str(pmid),
            'date': f"{date['Year']}-{date['Month']}-{date['Day']}",
            'journal': journal,
            'title': str(title),
            'abstract': abstract,
            'doi': str(doi),
            'authors': authors
        }

        paper_details.append(paper_info)

    print('abstract',len(paper_details), query_term)
    return paper_details

email = 'mutwil@gmail.com'
year_cutoff = 2005

papers = []

papers += search_and_fetch_journal_papers("Plant Cell[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plant Physiol[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plant J[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plant Cell Physiol[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plant Mol Biol[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Front Plant Sci[Journal]", email, 2020)
papers += search_and_fetch_journal_papers("J Exp Bot[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Mol Plant[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Planta[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plant Signal Behav[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("New Phytol[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("BMC Plant Biol[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Mol Plant Microbe Interact[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plant Cell Environ[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Nat Plants[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Physiol Plant[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("J Integr Plant Biol[Journal]", email, year_cutoff)
papers += search_and_fetch_journal_papers("J Plant Physiol[Journal]", email, year_cutoff)

papers += search_and_fetch_journal_papers("PNAS[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Science[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Cell[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("J Biol Chem[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("FEBS Lett[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Curr Biol[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("EMBO J[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Development[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Plos Genet[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Nature[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Genes Dev[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Nucleic Acid Res[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Biochem J[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Int J Mol Sci[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Nat Commun[Journal] AND plants[MeSH Terms]", email, year_cutoff)
papers += search_and_fetch_journal_papers("Elife[Journal] AND plants[MeSH Terms]", email, year_cutoff)
counter = 0

with open(path+'abstracts', 'wb') as file:
    pickle.dump(papers, file)

alreadyDone = os.listdir(path+'/annotations')

for i in papers:
    pubmed, abstract = str(i['pmid']), i['abstract']
    
    if pubmed+'.txt' not in alreadyDone:
        try:
            if abstract!='':
                result = gpt_analyzer(abstract)
                v = open(path+'/annotations/%s.txt' % pubmed,'w')
                v.write(str(abstract)+'\n\n'+str(result))
                v.close()
            else:
                v = open(path+'/annotations/%s.txt' % pubmed,'w')
                v.write('')
                v.close()
        except:
            print('gpt failed', i)
    else:
        print(pubmed, 'already done')
    counter +=1
    print(counter, 'out of', len(papers))
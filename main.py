__author__ = 'maximkuleshov'

import json
from itertools import combinations, product
from os import listdir
from os.path import isfile, join
from time import sleep

import requests


def pairs(*lists):
    for t in combinations(lists, 2):
        for pair in product(*t):
            if pair[0] != pair[1]:
                yield pair


def get_enrichr_results(gene_set_library, genelist, description):
    addlist_url = 'http://amp.pharm.mssm.edu/FlyEnrichr/addList'
    payload = {
        'list': (None, genelist),
        'description': (None, description)
    }

    response = requests.post(addlist_url, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    sleep(1)
    data = json.loads(response.text)

    query_string = 'http://amp.pharm.mssm.edu/FlyEnrichr/enrich?userListId={0}&backgroundType={1}'
    user_list_id = data['userListId']
    response = requests.get(query_string.format(user_list_id, gene_set_library))
    sleep(1)
    return [data['shortId'], json.loads(response.text)]


def get_libraries():
    # Returns names of all libraries
    libs_json = json.loads(requests.get('http://amp.pharm.mssm.edu/FlyEnrichr/datasetStatistics').text)
    libs = [lib['libraryName'] for lib in libs_json['statistics']]
    return libs


def parse_gmt(gmt):
    # Parses gmt-file to dictionary with format {term1: [gene1, gene2, ..., geneN], term2: [...], ...}
    gmt_dict = dict()
    for line in gmt:
        term, desc, *genes = line.strip().split('\t')
        genes = [gene.split(',')[0] for gene in genes]
        gmt_dict[term] = sorted(list(set(genes)))
    return gmt_dict


def main():
    # Assumption - gmt files names are the same as FlyEnrichr libraries names
    gmt_dir = 'gmt'
    # List of libraries based on content of 'gmt' folder
    gmt_libs = [f.split('.')[0] for f in listdir(gmt_dir) if isfile(join(gmt_dir, f))]
    # List of libraries form the FlyEnrichr website
    fly_libs = get_libraries()
    # Pairwise combinations of gmt libs and website libs. Skips pairs with duplicate libs
    gmt_fly_pairs = set(pairs(gmt_libs, fly_libs))

    for gmt_lib, fly_lib in gmt_fly_pairs:
        print('{0} vs {1}'.format(gmt_lib, fly_lib))
        gmt = parse_gmt(open('{0}/{1}.txt'.format(gmt_dir, gmt_lib), 'r').readlines())
        for term in gmt:
            print(term)
            genelist = '\n'.join(gmt[term])
            result = get_enrichr_results(fly_lib, genelist, '{0} against {1}'.format(term, fly_lib))
    return None


if __name__ == '__main__':
    main()

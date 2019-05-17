'''
Created on 23.02.2017

@author: marisakoe


If there are synonyms in the data, sample the words accordingly to take them into account.
This is done concept wise and only for languages which have synonyms in the data, the rest of the data will remain the same.
'''

import codecs, numpy

from collections import defaultdict

import helper_functions
from itertools import combinations, permutations, product



def sample_data(f):
    '''
    samples the data in case there are synonyms.
    If the list of words is longer than 1: sample
    Otherwise leace everything as it is
    :param f: input file
    '''
    ##data dict key=concept value=dict with key=lang value=list of words
    #data_ielexAsjp = helper_functions.read_IelexASJP(f)
    data = helper_functions.read_NlexIELex(f)
    new_data_dict = defaultdict(dict)
    ##for each concept in the dictionary, get the word list, get the length of the sublists, check if any value is greater than 1
    for con, lang_dict in data.items():
        if con == "Berg":
            ##the word list of this concept
            word_list_concept = lang_dict.values()
            #print word_list_concept
            ##a list with the lengths of the sublists to check for synonyms
            list_lenSublists = [len(x) for x in word_list_concept]
            #print list_lenSublists
            #matches = [x for x in list_lenSublists if x > 1]
            #print matches
            #number_samples = reduce(lambda x, y: x*y, matches)
            #print number_samples
            ##True if any character is greater than 1 otherwise False
            bol = any(i > 1 for i in list_lenSublists)
            if bol == True:
                ##sample according to the numbers of synonyms
                sample_dict = sample_wordlist(lang_dict)
                for iteration, dict_langs in sample_dict.items():
                    new_data_dict[con][str(iteration)] = dict_langs
      
            else:
                for lang, wrdList in lang_dict.items():
                    lang_dict[lang] = wrdList[0]
                new_data_dict[con][str("none")]=lang_dict
                #do nothing an include the concept as it is into the dictionary


    return new_data_dict


###########################helper methods#####################################
    
def sample_wordlist(lang_dict):
    '''
    sample the data depending on the product of the synonyms. Each combination should be present once.
    :param lang_dict: the dictionary of the languages in the concept
    :return new_dict: a dict with key=numer of iteration value=dict with key=lang value=single word in a list
    '''
    
    ##cerate an empty list for the samples
    sample_List = []
    ##for lang and wordlist in the dictionary
    for lang, wdList in lang_dict.items():
        ##check if length of wordlist > 1
        if len(wdList) > 1:
            print lang, wdList
            ##set a new list
            index_list = []
            ## get the index and item in the list
            for i,item in enumerate(wdList):
                ##append the language and the index of the word
                index_list.append((lang,i))
            ##append it to the bigger list, including list with language and index
            sample_List.append(index_list)
    ##set a new list for the combinations
    print sample_List
    combination_list = []
    ##build the product of the list of lists and append it to the list
    for tup in product(*sample_List):
        combination_list.append(tup)
    print len(combination_list)
    ##create a new dictionary for returning all dictionaries with possible combinations
    new_dict = defaultdict()

    ##for index, combination in the list
    for idx,com in enumerate(combination_list):
        ##create a dict for this combination
        l_dict = defaultdict()
        ## for each pair in the combination
        for pair in com:
            ##get the language
            lang = pair[0]
            ##get the index of the word
            idx_wd = pair[1]
            ##get the wordlist from the dictionary
            wdList = lang_dict[lang]
            ##append the language and the word to the dictionary
            l_dict[lang]=wdList[idx_wd]
        ##for all langs and wordlists in the dictionary
        for lang1, wdList1 in lang_dict.items():
            ##check if lang is not already in the dictionary and append the lang with only one word
            if not lang1 in l_dict:
                l_dict[lang1] = wdList1[0]
        ##fill the dictionary for each combination with the sample dictionary
        new_dict[idx+1] = l_dict
  
#     for x,y in new_dict.items():
#         print x,y
    return new_dict



if __name__ == '__main__':
    #f = codecs.open("../input/IELex-2016.tsv.asjp","r", encoding="utf-8")
    f = codecs.open("../input/ielex-northeuralex-intersection-with-asjp.tsv.asjp","r", encoding="utf-8")
    sample_data(f)
    f.close()
    
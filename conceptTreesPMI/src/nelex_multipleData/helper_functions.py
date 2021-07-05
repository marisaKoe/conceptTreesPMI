'''
Created on 31.07.2017

@author: marisa
'''

import codecs
from collections import defaultdict


def read_nelexPMI(data):
    '''
    reads the nelex datafile, which can be downloaded from the website with the ASJP transcriptions
    :param f:
    '''
    
    #initialize a default dict with the data
    d = defaultdict(lambda: defaultdict(list))
    cog_dict = defaultdict(lambda: defaultdict(list))
    data.readline()

    for line in data:
        line = line.strip()
        #split it by tab
        arr = line.split("\t")
        concept = arr[0]
        lang=arr[1]
        asjp=arr[2]
        cc = arr[-1]
        
        if " " in concept:
            concept = concept.replace(" ","")
        
        #word transcription
        asjp_word = arr[2].split(",")[0]
        asjp_word = asjp_word.replace(" ", "")
        #tokenized_word = ipa2tokens(asjp_word)
        #asjp_word = "".join(tokens2class(tokenized_word, 'asjp'))
        asjp_word = asjp_word.replace("%","")
        asjp_word = asjp_word.replace("~","")
        asjp_word = asjp_word.replace("*","")
        asjp_word = asjp_word.replace("$","")
        asjp_word = asjp_word.replace("\"","")
        asjp_word = asjp_word.replace(""" " ""","")
        asjp_word = asjp_word.replace("K","k")
        if len(asjp_word) < 1:
            continue
        if "K" in asjp_word:
            print line
            
#         ##check if the language in not in the dictionary yet, create the entry with a list of the first asjp_word
#         if not lang in d[concept]:
#             d[concept][lang] = [asjp_word]
#         ##if the language is already in the dictionary, append the list of asjp words
#         else:
        d[concept][lang].append(asjp_word)
        cog_dict[concept][lang].append((asjp_word,cc))

    

        
    return d, cog_dict


def writePhy(dm,f):
    '''
    Writes the distance matrix into a file of format phylip
    :param dm: the dictionary with the distances
    :param f: the output file name
    '''
    ntaxa = len(dm)
    f.write(str(ntaxa)+"\n")
        
    i=0
    for ka in dm:
        row=ka+" "
        for kb in dm:
            row=row+" "+str(dm[ka][kb])
        f.write(row+"\n")
    f.close()


if __name__ == '__main__':
    data = "../input/nelexAsjp.cognates"
    f = codecs.open(data,"r", encoding="utf-8")
    #returns one dictionary
    #key=concept value=dict with key=number iteration or non value=dict with key=lang value=list with one word (needs to be extracted)
    te = read_nelexPMI(f)
    f.close()
    
    
    
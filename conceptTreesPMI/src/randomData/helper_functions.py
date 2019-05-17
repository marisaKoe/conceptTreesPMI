'''
Created on 17.11.2016

@author: marisa
'''

import codecs
from collections import defaultdict

def read_data(f):
    '''
    reads the random data from a file
    :param f: filename
    :return d: a dictionary key = concept value = dict with key = lang+_+word value = word
    '''
    #initialize the default dictionary
    d = defaultdict(dict)
    f.readline()
    #set a counter for the concepts
    concept_counter = 0
    #go through the file
    for line in f:
        #increase the concept counter
        concept_counter += 1
        #split and strip the line
        line = line.strip().split("\t")
        #set the language counter
        language_counter = 0
        #go through the list of words
        for word in line:
            #increase the language counter
            language_counter += 1
            #set the dictionary
            d["c"+str(concept_counter)]["l"+str(language_counter)] = word
        
    
    return d


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
    f = codecs.open("../input/Random_words_10.txt","r", encoding="utf-8")
    read_data(f)
    f.close()
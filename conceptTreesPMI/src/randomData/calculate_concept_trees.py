'''
Created on 17.11.2016

@author: marisa
'''

from numpy import *
import itertools as it
from collections import defaultdict
import codecs, multiprocessing
from sklearn import metrics

from helper_functions import *


unique_chars = []
gp1 = -2.49302792222
gp2 = -1.70573165621

lodict={}

def nw(x,y,lodict=lodict,gp1=gp1,gp2=gp2):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    #length of the words
    n,m = len(x),len(y)
    dp = zeros((n+1,m+1))
    pointers = zeros((n+1,m+1),int)
    for i in xrange(1,n+1):
        dp[i,0] = dp[i-1,0]+(gp2 if i>1 else gp1)
        pointers[i,0]=1
    for j in xrange(1,m+1):
        dp[0,j] = dp[0,j-1]+(gp2 if j>1 else gp1)
        pointers[0,j]=2
    for i in xrange(1,n+1):
        for j in xrange(1,m+1):
            match = dp[i-1,j-1]+lodict[x[i-1],y[j-1]]
            insert = dp[i-1,j]+(gp2 if pointers[i-1,j]==1 else gp1)
            delet = dp[i,j-1]+(gp2 if pointers[i,j-1]==2 else gp1)
            dp[i,j] = max([match,insert,delet])
            pointers[i,j] = argmax([match,insert,delet])
    alg = []
    i,j = n,m
    while(i>0 or j>0):
        pt = pointers[i,j]
        if pt==0:
            i-=1
            j-=1
            alg = [[x[i],y[j]]]+alg
        if pt==1:
            i-=1
            alg = [[x[i],'-']]+alg
        if pt==2:
            j-=1
            alg = [['-',y[j]]]+alg
    return dp[-1,-1]


def compute_random_dm(data):
    '''
    main function to create distance matrices and write them into files
    :param data: the name of the data file
    '''
    #reads the file with the 41 sounds of ASJP
    f = open('input/sounds41.txt')
    sounds = array([x.strip() for x in f.readlines()])
    f.close()

    #reads the log odds scores for the whole word languages
    f = open('input/pmi-world.txt','r')
    l = f.readlines()
    f.close()
    logOdds = array([x.strip().split() for x in l],double)
    
    #assigns the log odds to the sound pairs to create the lodict
    for i in xrange(len(sounds)):#Initiate sound dictionary
        for j in xrange(len(sounds)):
            lodict[sounds[i],sounds[j]] = logOdds[i,j]

    

    concept_dict = defaultdict(lambda: defaultdict(int))
    #reads the data file
    f = codecs.open(data,"r", encoding="utf-8")
    #returns two dictionaries
    #te = key:concept and language value:word
    #cogid_dict = key=concept and language value:cognate class
    te = read_data(f)
    f.close()
    
    
    
    
    #for each concept in the dictionary
    for concept,langs in te.items():
        #print concept, langs
        dist_dict = defaultdict(lambda: defaultdict(float))
        langs_list = []
        #for each language in the dictionary, append it to a list
        for lang in langs:
            if lang not in langs_list:
                langs_list.append(lang)

        #for each language in the list, make pairs of languages
        #append the language pair to a list
        #append the language pair and their corresponding words to a list
        for lang_pair in it.combinations(langs_list,r=2):
            #print lang_pair
            l1, l2 = lang_pair
            #lang_pair_list.append((l1, l2))
            #directly computing pmi distance here
            wrd_score = 1.0-((2.0*nw(langs[l1],langs[l2]))/(nw(langs[l1],langs[l1])+nw(langs[l2],langs[l2])))
            #print nw(langs[l2],langs[l2])
            dist_dict[l1][l2] = wrd_score   
            dist_dict[l2][l1]=wrd_score 
        #for k,v in dist_dict.items():
        #    print k,v
        
        fout1 = codecs.open("output/randomData/phylip10/randomData_"+concept+".phy","wb","utf-8")
        writePhy(dist_dict, fout1)
        print "still creating matrices"
        #fout2 = codecs.open("output/nexus/"+out+"_"+concept+".nex","wb","utf-8")
        #writeNexus(dist_dict, fout2)

if __name__ == '__main__':
    data = "../input/Random_words_10.txt"
    compute_random_dm(data)
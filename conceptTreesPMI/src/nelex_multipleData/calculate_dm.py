'''
Created on 31.07.2017

@author: marisa

calculates the distance matrix with PMI (and sigmoid if needed).

For each concept:

+ compute dm
+ save mtx
+ if number of mtx > 1:
    + compute sdm
    + reconstruct tree
    + delete mtx
+ otherwise:
    + reconstruct single tree
'''


import codecs, multiprocessing, os
from numpy import *
import itertools as it
from collections import defaultdict
from sklearn import metrics


import reconstruct_trees
from helper_functions import read_nelexPMI,writePhy



######global variables#################
gp1 = -2.49302792222
gp2 = -1.70573165621

lodict={}




################helper methods######################

def compute_dm(data, method):
    '''
    main function to create distance matrices and write them into files
    :param data: the name of the data file
    '''
    print "in compute_dm(data)"
    #reads the file with the 41 sounds of ASJP
    f = open('input/sounds.txt')
    sounds = array([x.strip() for x in f.readlines()])
    f.close()
    
    #reads the log odds scores for the whole word languages
    f = open('input/pmi_matrix_nelex.txt','r')
    l = f.readlines()
    f.close()
    logOdds = array([x.strip().split() for x in l],double)
    
    #assigns the log odds to the sound pairs to create the lodict
    for i in xrange(len(sounds)):#Initiate sound dictionary
        for j in xrange(len(sounds)):
            lodict[sounds[i],sounds[j]] = logOdds[i,j]

    
    
    #reads the data file
    f = codecs.open(data,"r", encoding="utf-8")
    #returns one dictionary
    #key=concept value=dict with key=number iteration or non value=dict with key=lang value=list with one word (needs to be extracted)
    te, cogid_dict = read_nelexPMI(f)
    f.close()
    
    out = data.split(".")
    out = out[0].split("/")[-1]
    
    ############take synonyms into account by taking the wordpair with the smallest distance
    #for each concept in the dictionary
    for concept,langs_dict in te.items():           
                    
        dist_dict = defaultdict(lambda: defaultdict(float))
        langs_list = []
        #for each language in the dictionary, append it to a list
        for lang in langs_dict:
            if lang not in langs_list:
                langs_list.append(lang)

        #for each language in the list, make pairs of languages
        #append the language pair to a list
        #append the language pair and their corresponding words to a list
        for lang_pair in it.combinations(langs_list,r=2):
            l1, l2 = lang_pair
            word_list1 = langs_dict[l1]
            word_list2=langs_dict[l2]
            if len(word_list1)==1 and len(word_list2)==1:
                if method == "pmi":
                    ##directly computing pmi distance here
                    ##transform similarity value into distance value
                    wrd_score = 1.0-((2.0*nw(word_list1[0],word_list2[0]))/(nw(word_list1[0],word_list1[0])+nw(word_list2[0],word_list2[0])))
                
                elif method == "sigmoid":
                    ##sigmoid function zur Berechnung der Distanz
                    ##Vorteil Sigmoid ist eine S-Kurve. Alle Werte unter 0.5 signalisieren eine negative Aehnlichkeit (-> nicht aehnlich), all Werte ueber 0.5 
                    ## signalisieren eine positive Aehnlichkeit (-> aehnlich). Dadurch kann der threshold fuer andere funktionen (lingpy) immer 0.5 sein.
                    wrd_score = 1.0 - (1.0/(1 + exp(-nw(word_list1[0],word_list2[0]))))

                dist_dict[l1][l2] = wrd_score   
                dist_dict[l2][l1]=wrd_score 
            else:
                combination_list = list(it.product(word_list1, word_list2))
                list_score=[]
                for word_pair in combination_list:
                    if method == "pmi":
                        ##directly computing pmi distance here
                        ##transform similarity value into distance value
                        wrd_score = 1.0-((2.0*nw(word_pair[0],word_pair[1]))/(nw(word_pair[0],word_pair[0])+nw(word_pair[1],word_pair[1])))
                        list_score.append(wrd_score)
                        dist_dict[l1][l2] = min(list_score)   
                        dist_dict[l2][l1]= min(list_score) 
                        
                    elif method == "sigmoid":
                        ##sigmoid function zur Berechnung der Distanz
                        ##Vorteil Sigmoid ist eine S-Kurve. Alle Werte unter 0.5 signalisieren eine negative Aehnlichkeit (-> nicht aehnlich), all Werte ueber 0.5 
                        ## signalisieren eine positive Aehnlichkeit (-> aehnlich). Dadurch kann der threshold fuer andere funktionen (lingpy) immer 0.5 sein.
                        wrd_score = 1.0 - (1.0/(1 + exp(-nw(word_pair[0],word_pair[1]))))
                        list_score.append(wrd_score)
                        dist_dict[l1][l2] = min(list_score)   
                        dist_dict[l2][l1]= min(list_score) 
                
                
                
                
        fout1 = codecs.open("output/nelex_multipleData/"+method+"/phylip/"+out+"_"+concept+".phy","wb","utf-8")
        writePhy(dist_dict, fout1)
        print "still creating matrices"
                
                




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

if __name__ == '__main__':
    data = "../input/nelexAsjp.cognates"
    compute_dm(data)
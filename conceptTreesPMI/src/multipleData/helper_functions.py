'''
Created on 23.02.2017

@author: marisakoe



helper functions:

* read the data of IELex+ASJP
* read the data of NothEuraLex+IELex
* write matrices into files for other methods
'''
import codecs
from collections import defaultdict

unique_chars = []

def read_IelexASJP(f):
    '''
    read the ielex+asjp data and save it in a dictionary.
    All synonyms are taken into account.
    :param f: the input file
    :return d: the dictionary with all words key=concept value=dict with key=lang and value=list with all words
    :return cogid_dict: the dictionary with the corresponding cognate classes key=concept value=dict with key=lang and value=list with all cognates
    '''
    #initialize a default dict
    d = defaultdict(dict)
    cogid_dict = defaultdict(dict)
    f.readline()
    data=defaultdict(dict)
    unique_langs = defaultdict()
    
    #for each line in the file
    for line in f:
        line = line.strip()
        #split it by tab
        arr = line.split("\t")
        #language
        lang = arr[0]
        #iso-code
        iso = arr[1]
        #gloss
        gloss = arr[2]
        #concept number
        concept = arr[3]
        #local-id
        local_id=arr[4]
        #transcription
        trans = arr[5]
        #cognate class
        cogid = arr[6]
        #notes
        note=arr[-1]
        cogid = cogid.replace("-","")
        cogid = cogid.replace("?","")
        #word transcription
        asjp_word = arr[5].split(",")[0]
        asjp_word = asjp_word.replace(" ", "")
        #tokenized_word = ipa2tokens(asjp_word)
        #asjp_word = "".join(tokens2class(tokenized_word, 'asjp'))
        asjp_word = asjp_word.replace("%","")
        asjp_word = asjp_word.replace("~","")
        asjp_word = asjp_word.replace("*","")
        asjp_word = asjp_word.replace("$","")
        asjp_word = asjp_word.replace("\"","")
        asjp_word = asjp_word.replace("K","k")
        if len(asjp_word) < 1:
            continue
        if "K" in asjp_word:
            print line
        for x in asjp_word:
            if x not in unique_chars:
                unique_chars.append(x)
        ##check if the language in not in the dictionary yet, create the entry wiht a list of the first asjp_word
        if not lang in d[concept]:
            d[concept][lang] = [asjp_word]
        ##if the language is already in the dictionary, append the list of asjp words
        else:
            d[concept][lang].append(asjp_word)
        
        ##do the same for the cognates (the lists are ordered so the cogid is always on the same position than the word in the other dictionary.
        if not lang in cogid_dict[concept]:
            cogid_dict[concept][lang] = [(asjp_word,cogid)]
        else:
            cogid_dict[concept][lang].append((asjp_word,cogid))
            
#     for con, l in d.items():
#         if con == "639":
#             print l
#     
#     for c,l in cogid_dict.items():
#         if c=="639":
#             print l
    #return (d, cogid_dict)  
    return d  

def read_NlexIELex(f):
    '''
    read the NorthEuraLex+IELex data and create a dictionary.
    Take account of all the synonyms.
    :param f: the input file
    '''
    #initialize a default dict
    d = defaultdict(dict)
    cogid_dict = defaultdict(dict)
    f.readline()
     
    unique_langs = []
    #for each line in the file
    for line in f:
        line = line.strip()
        #split it by tab
        arr = line.split("\t")
        #print arr
        #iso-code (language name)
        iso = arr[0]
        #gloss global
        #gloss_global = arr[1]
        #gloss ielex
        #gloss_ielex = arr[2]
        #gloss northeuralex
        gloss_nlex = arr[3]
        gloss_nlex =  gloss_nlex.split(":")
        gloss_nlex = gloss_nlex[0]
        
        #ortho ielex
        #ortho_ielex=arr[4]
        #ortho northeuralex
        #ortho_nlex = arr[5]
        #phono northeuralex
        #phon_nlex = arr[6]
        #transcription
        asjp_word = arr[7]
        #print asjp_word
        #print asjp_word
        #cognate class
        cogid = arr[8]
        #loan judgement
        #loan=arr[-1]
 
        for x in asjp_word:
            if x not in unique_chars:
                unique_chars.append(x)
         
        if iso not in unique_langs:
            unique_langs.append(iso)
        
        
        #print gloss_nlex
        ## if the iso code is not in the dictionary, set the dictionary with a list of words
        if not iso in d[gloss_nlex]:
            d[gloss_nlex][iso] =  [asjp_word]
        ##otherwise append the words to the list
        else:
            d[gloss_nlex][iso].append(asjp_word)
             
        #cogid_dict[gloss_nlex][iso] = [asjp_word,cogid]
         
        ## if the iso coade is not in the dictionary, set the dictionary with a list of cogids and words
        if not iso in cogid_dict[gloss_nlex]:
            cogid_dict[gloss_nlex][iso] = [(asjp_word,cogid)]
        else:
            cogid_dict[gloss_nlex][iso].append((asjp_word,cogid))
            

              
    #return (d, cogid_dict)
#     for x,y in d.items():
#         print x,y
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
    #f = codecs.open("../input/IELex-2016.tsv.asjp","r", encoding="utf-8")
    #read_IelexASJP(f)
    f = codecs.open("../input/ielex-northeuralex-intersection-with-asjp.tsv.asjp","r", encoding="utf-8")
    read_NlexIELex(f)
    f.close()
    
    
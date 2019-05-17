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
    d = defaultdict(dict)
    ##cogante dicht
    cogid_dict = defaultdict(dict)
    
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
            #print concept
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
            
        #fill the dictionary with key=concept val=newDict with key=language+word val=word
        d[concept][lang] = asjp
        
        #fill the dictionary with key = concept and language and value = word transcription
        #d[concept+":"+lang] = asjp_word
        #fill the dictionary with key = concept and language and value = cognate identity
        cogid_dict[concept][lang] = cc
    

        
    return (d, cogid_dict)




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
    
    
def write_cognates(out, cogid_dict):
    '''
    writes the languages and cognate ids into a file for each concept
    :param out: the name ouf the inputfile
    :param cogid_dict:dictionary key=concept value=dict with key=language and value=cognate id
    '''
    for concept,langs in cogid_dict.items():
        fout1 = codecs.open("output/nelex/cognates/"+out+"_"+concept+".csv","wb","utf-8")
        fout1.write(concept+"\n")
        for lang, cogid in langs.items():
            fout1.write(lang+"\t"+cogid+"\n")
        fout1.close()

if __name__ == '__main__':
    data = "../input/nelexAsjp.cognates"
    f = codecs.open(data,"r", encoding="utf-8")
    #returns two dictionaries
    #te = key:concept and language value:word
    #cogid_dict = key=concept and language value:cognate class
    te, cogid_dict = read_nelexPMI(f)
    f.close()




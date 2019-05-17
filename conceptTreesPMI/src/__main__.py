'''
Created on 14.07.2016

@author: marisa
'''

#from ielex import main, reconstruct_trees_phy

#from IELexNLex import compute_dm, reconstruct_trees_phy

#from randomData import compute_random_dm, reconstruct_trees_phy_random

#from multipleData import compute_dm

#from nelex import compute_dm, reconstruct_trees_phy

from nelex_multipleData import compute_dm, reconstruct_trees_phy

if __name__ == '__main__':
    ####ielex and asjp from Taraka
    #data = "input/IELex-2016.tsv.asjp"
    #main(data)
    #reconstruct_trees_phy()
    
    
    ###ielex/NorthEuraLex and asjp from Johannes
    #data = "input/ielex-northeuralex-intersection-with-asjp.tsv.asjp"
    #compute_dm(data)
    #reconstruct_trees_phy()
    
    #######random Data ASJP
    #data = "input/Random_words_10.txt"
    #compute_random_dm(data)
    #reconstruct_trees_phy_random()
    
    ###multiple data including Synonyms 
    #data = "input/ielex-northeuralex-intersection-with-asjp.tsv.asjp"
    #compute_dm(data)
    
    ######nelex without synonyms
    #data = "input/nelexAsjp.cognates"
    #compute_dm(data)
    #reconstruct_trees_phy()
    
    ###nelex multiple data
    data = "input/nelexAsjp.cognates"
    methods=["pmi","sigmoid"]
    
    for method in methods:
        compute_dm(data, method)
        reconstruct_trees_phy(method)
    
    
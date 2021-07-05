'''
Created on 14.07.2016

@author: marisa
'''


#from nelex import compute_dm, reconstruct_trees_phy

from nelex_multipleData import compute_dm, reconstruct_trees_phy

if __name__ == '__main__':

    
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
    
    
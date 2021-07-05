'''
Created on 31.07.2017

@author: marisa

'''


import glob
import rpy2.robjects as r
from rpy2.robjects.packages import importr

from numpy import *
from helper_functions import writePhy

#import the basis of R
base = importr("base")
utils = importr("utils")
stats = importr("stats")
#imports ape, required for phangorn
ape = importr("ape")
#imports phangorn
phangorn = importr("phangorn")



def reconstruct_trees_phy(method):
    '''
    Reconstructs Neighbor Joining trees for each concept
    The trees are reconstructed using R and the packages ape and phangorn
    '''
    list_matrices = glob.glob("output/nelex_multipleData/"+method+"/phylip/*.phy")
    for f in list_matrices:
        out = f.split(".")
        out = out[0].split("/")
        out = out[-1]
        t = utils.read_table(f, skip=1, row_names=1)
        mx = base.as_matrix(t)
        dm = stats.as_dist(mx)
        #nj trees
        treeNJ = ape.nj(dm)
        ape.write_tree(treeNJ, file="output/nelex_multipleData/"+method+"/phylipTrees/nj/"+out+"+njTree.nwk")
        #fastme trees
        treeFastme = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
        ape.write_tree(treeFastme, file="output/nelex_multipleData/"+method+"/phylipTrees/fastme/"+out+"+fastmeTree.nwk")


        
        




if __name__ == '__main__':
    pass
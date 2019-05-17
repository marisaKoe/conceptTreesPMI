'''
Created on 24.02.2017

@author: marisakoe

1. Method for normal tree reconstruction


2. Method for SDM and consensus tree reconstruction

One method to create a consensus tree is the SDM method.
Normally, this method is used, if the distance matrices contain a different number of taxa.
The method can be used here too, although all distance matrices have the same number of taxa.
The reason:
+ a matrix is needed for the Goodman-Kruskal-jackknife method

Try also: normal consensus tree which can be used in jackknifeLoanwordDistance

Look for an alternative method.


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

def reconstruct_trees_phy(f, out):
    '''
    Reconstructs Neighbor Joining trees for each concept
    The trees are reconstructed using R and the packages ape and phangorn
    '''
    
    
    t = utils.read_table(f, skip=1, row_names=1)
    mx = base.as_matrix(t)
    dm = stats.as_dist(mx)
    #nj trees
    treeNJ = ape.nj(dm)
    ape.write_tree(treeNJ, file="output/multipleData/phylipTrees/nj/"+out+"+njTree.nwk")
    #fastme trees
    treeFastme = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
    ape.write_tree(treeFastme, file="output/multipleData/phylipTrees/fastme/"+out+"+fastmeTree.nwk")
        
        
def create_concept_consensus(matrices, seqLen, concept):
    '''
    create concept and consensus trees for one concept
    :param consensusTreeDict:a list for the consensus tree
    :param conceptTreeDict:a list for the concept trees
    :param matrices:a list with all matrix files
    :param seqLen:the number of matrices for the concept
    :return conceptTrees: a list with all newick strings for this concept
    :return consensusTree: a list with the consensus tree for this concept
    '''
    #print concept
    #concept is used for the naming of the output file
    #seqLen-1 is needed for the sdm function (same for each matrix)
    #matrices are the files to be read in
    #create a list with the indices needed for the SDM method (n-1)
    seqList = [0]*seqLen
    for idx,num in enumerate(seqList):
        seqList[idx] = seqLen
    #print "in create consensus"
    #create a list for all distance matrices as R objects
    all_dm = []
    #for each matrix in the list
    for mtx in matrices:
        #read it in
        t = utils.read_table(mtx,header=False, skip=1, row_names=1)
        
        #make a matrix
        mx = base.as_matrix(t)
        #make a distance matrix
        dm = stats.as_dist(mx)

        #append it to list
        all_dm.append(dm)
    num_rows = t.nrow
    #print type(NA)
    #combine both lists
    overall_list = all_dm + seqList
    #give the elements of a list into the SDM method from ape
    dd = ape.SDM(*overall_list)
    #print dd[0]
    #create a numpy array out of it
    #matrixArray = squeeze(asarray(dd[0]))
    #print matrixArray
    #write the matrix into a file
    #phylip_output(matrixArray, "output/multipleData/phylip/" + concept + ".phy")
    #write the consensus matrix in a file
    base.cat(str(num_rows)+"\n", file="output/multipleData/phylip/" + concept + ".phy")
    utils.write_table(dd[0], file = "output/multipleData/phylip/" + concept + ".phy",quote = False,row_names=True, col_names=False, append=True)
    #use the first matrix to create a tree (also try fastme)
    treeNJ = ape.nj(dd[0])
    #write the tree in a file
    ape.write_tree(treeNJ, file= "output/multipleData/phylipTrees/nj/"+concept+"+njTree.nwk")
    ##fastme trees
    treeFastme = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
    ape.write_tree(treeFastme, file="output/multipleData/phylipTrees/fastme/"+concept+"+fastmeTree.nwk")



if __name__ == '__main__':
    
    pass
    
    
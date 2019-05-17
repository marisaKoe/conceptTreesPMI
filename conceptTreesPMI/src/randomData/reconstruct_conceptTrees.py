'''
Created on 17.11.2016

@author: marisa
'''
import glob
import rpy2.robjects as r
from rpy2.robjects.packages import importr

#import the basis of R
base = importr("base")
utils = importr("utils")
stats = importr("stats")
#imports ape, required for phangorn
ape = importr("ape")
#imports phangorn
phangorn = importr("phangorn")

def reconstruct_trees_phy_random():
    '''
    Reconstructs Neighbor Joining trees for each concept
    The trees are reconstructed using R and the packages ape and phangorn
    '''
    #print "in reconstruct trees"
    list_matrices = glob.glob("output/randomData/phylip10/*.phy")
    #print list_matrices
    for f in list_matrices:
        out = f.split(".")
        #print out
        out = out[0].split("/")
        out = out[-1]
        t = utils.read_table(f, skip=1, row_names=1)
        mx = base.as_matrix(t)
        dm = stats.as_dist(mx)
        #nj trees
        tree = ape.nj(dm)
        ape.write_tree(tree, file="output/randomData/phylipTrees10/nj/"+out+"+njTree.nwk")
        #fastme trees
        tree1 = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
        ape.write_tree(tree1, file="output/randomData/phylipTrees10/fastme/"+out+"+fastmeTree.nwk")
        
        

if __name__ == '__main__':
    reconstruct_trees_phy_random()
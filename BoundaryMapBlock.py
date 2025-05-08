import numpy as np
import itertools
from sage.all import *
from dataclasses import dataclass
import multiprocessing
from functools import partial
from SignedPerm import *
from basis import *
def generateSignBlocks(n,Basis=None):
    if Basis==None:
        Basis=matrix(matrix.identity(2**(n-1)),base_ring =GF(Integer(2)))
    print(Basis)
    N=list(range(0,n))
    trans=itertools.combinations(N,2)
    diags=generateDiagTurOr(n)
    transitionBlocks=dict()
    for tran in trans:
        transign0=[1 for _ in range(0,n)]
        transign1=[1 for _ in range(0,n)]
        transign0[tran[0]]=-1
        transign1[tran[1]]=-1
        matrixDict=dict()
        for i in range(len(diags[1])):
            for j in range(len(diags[-1])):
                mult= signMul(diags[1][i],diags[-1][j])
                matrixDict[(i,j)]=int( mult==transign0 or mult == transign1  )
        transitionBlocks[tran] = Basis.inverse()*matrix(matrixDict,base_ring=GF(Integer(2)))*Basis

    return transitionBlocks
def getBndryMapBlock(domB,ranB,n):# dom and range are permutation not signed permutations 
    matDict=dict()
    for d in domB:
        for r in ranB:

            rinv=permInv(r)
            res=permMul(d,rinv)
            if is_tran(res):
                matDict[ (r,d)]= permTotran(res)
            else:
                matDict[ (r,d)]= 0 
    return matDict
def getBndryMapBlocks( n):
   perms= generatePermutationsInOrder(n,permRLex)
   matdicts=[]
   for i in range(1,len(perms)):
       matdicts.append(  (getBndryMapBlock(perms[i-1],perms[i],n),perms[i-1],perms[i] ) )
   return matdicts
def getBndryMapByBlocks(n,Basis=None):
    #perms= generatePermutationsInOrder(n,permRLex)
    blocks=getBndryMapBlocks(n)
    signBlocks=generateSignBlocks(n,Basis=Basis)
    bndry=dict()
    for i in range(0,len(blocks)):
        print(i)
        domB=blocks[i][1]
        ranB=blocks[i][2]
        blockmatrix=[[signBlocks[blocks[i][0][(r,d)]] for d in domB ] for r in ranB ]
        bndry[i+1] = block_matrix(blockmatrix)
    return bndry




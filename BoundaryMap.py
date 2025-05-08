import numpy as np
import itertools
from sage.all import *
from dataclasses import dataclass
import multiprocessing
from functools import partial
from SignedPerm import *
from basis import *

def getindexOfOne(ranB,n,Trans,domB):

    listofone=[]

    for i in range(len(domB)):

        for j in range(len(ranB)):
            for Tran in Trans: 
                res=signPerm.mul(Tran,domB[i])
                if signPerm.is_equal(res,ranB[j]):
                    listofone.append((i,j))
                    break
    return listofone

def getBndryMap(domB,ranB,Trans,n):
    pool = multiprocessing.Pool(16)
    getindexOfOneP=partial( getindexOfOne,ranB,n,Trans)
    listofone = pool.imap_unordered(getindexOfOneP,(domB,))
    boundarymap= [[  0  for i in range(len(domB))] for j in range(len(ranB))]
    for pairs in listofone:
        for pair in pairs:
            i=pair[0]
            j=pair[1]
            boundarymap[j][i]=1
    return matrix(boundarymap,base_ring=GF(Integer(2)))
def getBndryMaps(n,gB=generateBasis):
    Trans=generateSignedTrans(n)
    basis=gB(n)
    bndryMaps=dict()
    for i in range(1, n*(n-1)//2+1):
        print(str(i) + "th boundarymap is being calculated")
        bndryMap=getBndryMap(basis[i-1],basis[i], Trans,n)
        bndryMaps[i]=(bndryMap)
    return bndryMaps

import numpy as np
import itertools
from sage.all import *
from dataclasses import dataclass
def permMul(p1,p2): # multiplies to permutations given as a list reprensenting oneline notations
    res=[]
    for i in range(len(p1)):
        res.append(p1[p2[i]])
    return res
def permInv(p):  #O(n^2) metod for inverting a permutation given as a list reprensenting oneline notation
    n=len(p)
    pinv=[]
    for i in range(n):
        pinv.append(p.index(i))
    return pinv
def permAct(perm,sign): # action of a permutation to a list. As a matrix this is conjugating a sign matrix by permutation
    res=sign.copy()
    for i in range(len(sign)):
        res[i]=sign[perm[i]]
    return res
def signMul(s1,s2): # multiplies 2 sign array. 
    res=[]
    for i in range(len(s1)):
        res.append(s1[i]*s2[i])
    return res
def diagMat(diagVal): # returns a diagonal matrix with given diagonal entries, takes a length n array and returns n by n matrix
    n=len(diagVal)
    Arr= [ [ 0 for i in range(0,n)]  for j in range(0,n) ]
    for i in range(0,n):
        Arr[i][i]=diagVal[i]

    return np.asarray(Arr)
def arrkthEntry(dim,k): # returns a lengtn dim array with k th entry 1 and rest zero
    A=dim*[0]
    A[k]=1
    return A
def permTotran(p):# given a transposition in oneline notation as array turns it to trasnposition notation with tuples 
    for i in range(len(p)):
        if i!=p[i]:
            return (i,p[i])
def is_tran(p): # returns true if given permutations is a transposition
    n=len(p)
    count=0
    for i in range(n):
        if i!=p[i]:
            count+=1
    if count == 2:
        return True
    else:
        return False

def permTomat(permArr,dim): # returns permutation matrix for given permutations (array reprensenting oneline notation)

    return np.asarray( [ arrkthEntry(dim,a) for a in permArr  ])
@dataclass(frozen=True)
class signPerm: # class for signed permutation perm and its sign.
    perm:list
    sign:list
    size:int
#    def __init__(self,perm,sign,n):
#        self.size=n
#        self.perm=perm
#        self.sign=sign
    def mul(sp1,sp2): # sign permutation is considered to be sign matrix times permutation matrix. This implement multiplications of two signed permutations using  comutation relations. Calculate result without using matrix
        perm1=sp1.perm.copy()
        perm2=sp2.perm.copy()
        sign1=sp1.sign.copy()
        sign2=sp2.sign.copy()
        resperm= permMul(perm1,perm2)
        ressign= signMul(sign1,permAct(perm1,sign2))
        res=signPerm(resperm,ressign,sp1.size)
        return res
    def is_equal(sp1,sp2):# return true if two given matrix are equal.

        perm1=sp1.perm.copy()
        perm2=sp2.perm.copy()
        sign1=sp1.sign.copy()
        sign2=sp2.sign.copy()
        return (perm1==perm2) and (sign1==sign2)
    def matrix(self): # returns matrix represatation of the sign matrix
        permMat=permTomat(self.perm,self.size)
        signMat=diagMat(self.sign)
        res=np.matmul(signMat,permMat)
        return matrix(res)
    def copy(self): # not needed return itself as new objeck
        return signPerm(self.perm,self.sign,self.size)

def countInv(perm): # counting inversion O(n^2) can be improved to O(nlogn) or python built-in functions maybe implemented 
    inv=0
    for i in range(0,len(perm)):
        for j in range(i+1, len(perm)):
            if perm[j]<perm[i]:
                   inv+=1 # inv=inv+1
    return inv
def permLex(perm): # assign a number to a permutation to give lexicographical order. Used as a order key to make lexicographical
    strord=''
    for p in perm:
        strord+=str(p)
    return int(strord)
def permRLex(perm):# reverse order of above
    strord=''
    for p in perm:
        strord+=str(p)
    return -int(strord)


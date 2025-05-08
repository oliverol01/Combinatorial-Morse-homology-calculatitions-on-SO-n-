from SignedPerm import *
import itertools

def generatePermutationsInOrder(n,orderKey=None, as_list=False):
    permList= [[] for _ in range(0,n*(n-1)//2 +1)]
    for perm in itertools.permutations(range(0,n)):
        if as_list==True:
            perm=list(perm)
        temp=permList[countInv(perm)]
        temp.append(perm)
    if orderKey==None:
        return permList
    permListSorted=[]
    for i in range(0,n*(n-1)//2 +1):
        tmp=sorted(permList[i],key=orderKey)
        permListSorted.append(tmp)
    return permListSorted
def permReorderbytop(permList,n):# reorder the   permutations with inversion more than L/2= n*(n-1)/4 using bijective map  multiplication by (n,n-1,......,1) from left between k and L-k s 
    L=len(permList)-1
    
    mid=L//2
    top= list(reversed(range(0,n)))
    for k in range(0,mid+1):
        permList[L-k]= [permMul( perm,top)  for perm in permList[k]]
    return permList
def generateDiagTurOr(n,k=0,Last=None):
    if Last==None:
        LastC=[[1]]
    else:
        LastC=Last.copy()
    j=k
    if k==n-1:
        LastMin=deepcopy(LastC)

        for i in range(len(LastC)):
            LastMin[i][1]=-LastMin[i][1]   
        return deepcopy(dict({-1:LastMin, 1:LastC,}))
    else:
        for i in range(len(LastC)):
            LastC[i].append(1)
        for i in range(len(LastC)):
            tmp=LastC[i].copy()
            tmp[k],tmp[k+1]=  -tmp[k],-tmp[k+1]
            LastC.append(tmp)
        return generateDiagTurOr(n,j+1,LastC)
def halfDiagsTuOr(diags):
    minusdiag=[]
    plusdiag=[]
    for i in range(0,len(diags[-1]),2):
        minusdiag.append(diags[-1][i])
    for i in range(0,len(diags[1]),2):
        plusdiag.append(diags[1][i])
    return deepcopy(dict({-1:minusdiag, 1:plusdiag}))
def generateBasis(n,reorderByTop=False,half=False):
    perms=generatePermutationsInOrder(n,permRLex,as_list=True)
    if reorderByTop==True:
        perms=permReorderbytop(perms,n)
    diags=generateDiagTurOr(n)
    if half==True:
        diags=halfDiagsTuOr(diags)
    basis=[[] for _ in range(n*(n-1)//2 + 1) ]
    for inv in range(0,n*(n-1)//2 + 1):
        for perm in perms[inv]:
            for diag in diags[(-1)**inv]:
                basis[inv].append(signPerm(perm,diag,n))
    return basis

        
def generateBasisErol(n):
    perms=generatePermutationsInOrder(n,permRLex,as_list=True)
    diags=generateDiagTurOr(n)
    basis=[[] for _ in range(n*(n-1)//2 + 1) ]
    
    for inv in range(0,n*(n-1)//2 + 1):
        for perm in perms[inv]:
            permtoact=perm
            if (-1)**inv == 1:
                permtoact=permInv(perm)
            diagsreordered= [ permAct(permtoact ,diag) for diag in  diags[(-1)**inv]]
            print(diagsreordered)
            for diag in diagsreordered :
                basis[inv].append(signPerm(perm,diag,n))
    return basis

    
def generateTrans(n):
    N=range(0,n)

    trpsOL= []
    trps=itertools.combinations(N, 2)

    for trp in trps:
        trpOL=range(0,n)
        i=trp[0]
        j=trp[1]
        trpOL[i],trpOL[j] = trpOL[j],trpOL[i]
        trpsOL.append(trpOL)

    return trpsOL
def generateSignedTrans(n):
    N=range(0,n)

    trpsOL= []
    trps=itertools.combinations(N, 2)

    for trp in trps:
        trpOL=[i for i in range(0,n)]
        i=trp[0]
        j=trp[1]
        sign1=[2*int(k!=i) -1  for k in range(0,n)]
        sign2=[ 2*int(k!=j) -1 for k in range(0,n)]
        trpOL[i],trpOL[j] = trpOL[j],trpOL[i]
        #print(trpOL)
        signedTrp1=signPerm(trpOL,sign1,n)
        signedTrp2=signPerm(trpOL,sign2,n)
        trpsOL.append(signedTrp1)
        trpsOL.append(signedTrp2)
    return trpsOL
def generateTuranBasisM(n):
    first=2**{n-1}*[0]
    first[0]=1



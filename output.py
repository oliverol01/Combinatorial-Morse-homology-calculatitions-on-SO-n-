from SignedPerm import *
from basis import *
from BoundaryMap import *
from BoundaryMapBlock import *
import numpy as np
import itertools
from sage.all import *
import subprocess


def printMatBlock(matDict,domB,ranB,html_file):
    html = "<table border='1' style='border-collapse: collapse; text-align: center;'>\n"

# Add header row (domB)
    html += "  <tr><th></th>" + "".join(f"<th>{col}</th>" for col in domB) + "</tr>\n"

# Add table rows (ranB)
    for row in ranB:
        html += f"  <tr><th>{row}</th>"  # ranB header
        for col in domB:
            value = matDict.get((row, col), "")

            html += f"<td>{value}</td>"
        html += "</tr>\n"

# Close the table
    html += "</table>"

# Save to an HTML file
    with open(html_file, "w") as file:
        file.write(f"<html><body>{html}</body></html>")
def printSo(n,showmatrix=False,Basis=None):
    latexfile="\\documentclass[8pt]{article}"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\newcommand{\\Z}{\\mathbb{Z}}\n"
    latexfile+= "\\usepackage[a2paper,"
    latexfile+=        "bindingoffset=0.2in,"
    latexfile+=         "left=0in,"
    latexfile+=        "right=0.5in,"
    latexfile+=        "top=1in,"
    latexfile+=        "bottom=1in,"
    latexfile+=        "footskip=.25in]{geometry}\n"
    latexfile+="\\begin{document}\n"

    latexfile+=(f" \\( So( {n} )\\)\n")
    
    bmsage=getBndryMapByBlocks(n,Basis=Basis) 
    print(bmsage.keys())
    for l in range(1,n*(n-1)//2+1 ):
        L=n*(n-1)//2+1
        pair=L-l
        noteq=""
        if bmsage[pair] != bmsage[l].transpose() :
            noteq="not"
        latexfile +=  f" {l}th boundary map is "+noteq+   f" transpose of {pair}\n"
        latexfile+= "\\newline\n"


    C=ChainComplex(bmsage,base_ring=GF(Integer(2)))
    Hom= C.homology(generators=True)
    Betti=C.betti()
    for k in Hom.keys():
        tmpHom=Hom[k]
        latexfile+=f"{k} th homology generators are \n"
        for gen in Hom[k]:
            latexfile += "\\[" +  latex(gen[1].vector(k)) + "\\]\n"
    for k in Betti.keys():
        latexfile+=f"\\[H_{k} = \\Z_2^{Betti[k]} \\]\n"
    if showmatrix==True:
        for l in range(len(bmsage)):
            latexfile +=  f" {l}th boundary map is"
            latexfile+= "\\[ " + latex(bmsage[l] )  + "\\]\n"

    latexfile+="\\end{document}"
    file_name = f"So{n}.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 

def printSignTransBlock(n,Basis=None):
    sb=generateSignBlocks(n,Basis=Basis)
    N=range(0,n)
    trans=itertools.combinations(N,2)

    latexfile="\\documentclass[8pt]{article}"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\begin{document}\n"
    for t in trans:
        
        latexfile+="\\[ \\sigma_{" + f"{t}"+"}=" +  latex(sb[t])  + "  \\]\n"
        if sb[t].transpose()==sb[t]:
            latexfile+=" \\(\\sigma_{" + f"{t}"+"} \\) is symmetric"
    latexfile+="\\end{document}"
    file_name = f"So{n}SignBlocks.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 
def printKernelOfSignTransBlock(n):
    sb=generateSignBlocks(n)
    latexfile="\\documentclass[8pt]{article}"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\begin{document}\n"
    for k in sb.keys():
        latexfile+="\\( \\sigma_{" + f"{k}"+"}\\)"+"'s kernel is span of   \n"
        latexfile+="\\begin{align*}"
        a = sb[k].right_kernel().basis_matrix()
        totlen=0
        for b in a.rows():
            totlen+= len(b)
            latexfile+=latex(b) + ","
            if totlen ==32:
                latexfile+= "\\\\ \n"
                totlen=0
 
        latexfile = latexfile[:-1]
        latexfile += "\\end{align*}\n"
    latexfile+="\\end{document}"
    file_name =f"So{n}SignBlocksKernel.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 

def printImageOfSignTransBlock(n):
    sb=generateSignBlocks(n)
    latexfile="\\documentclass[8pt]{article}"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\begin{document}\n"
    for k in sb.keys():
        latexfile+="\\( \\sigma_{" + f"{k}"+"}\\)"+"'s image is span of   \n"
        latexfile+="\\begin{align*}"
        a = sb[k].column_space().basis()
        totlen=0
        for b in a:
            totlen+= len(b)
            latexfile+=latex(b) + ","
            if totlen ==32:
                latexfile+= "\\\\ \n"
                totlen=0
 
        latexfile = latexfile[:-1]
        latexfile += "\\end{align*}\n"
    latexfile+="\\end{document}"
    file_name=f"So{n}SignBlocksImage.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 

def listAllzeroSumOfblocks(n):
    sb=generateSignBlocks(n)
    for r in range(2,len(sb.keys())+1):
        for  ks in itertools.combinations(sb.keys(),r):
            s=sum([sb[k] for k in ks])
            if s==0:
                print(ks)
def printImageOfBndry(n):
    bm=getBndryMaps(n)
    latexfile="\\documentclass[8pt]{article}\n"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\begin{document}\n"

    latexfile+= "\\centering\n"

    for k in bm.keys():
        latexfile+=f"{k}th Boundary Map's image is spanned by\n "
        latexfile+="\\begin{align*}"
        count=0

        for b in bm[k].column_space().basis():
            latexfile+=latex(b.column()) +",\n"
            count+=1
            if count==10:
                latexfile += "\\end{align*}\n"
                latexfile+="\\begin{align*}\n"
                latexfile+="\\newline"

                count=0
                
        latexfile += "\\end{align*}\n"
    latexfile+="\\end{document}"
    file_name=f"So{n}BndryImage.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 
def printKernelOfBndry(n):
    bm=getBndryMaps(n)
    latexfile="\\documentclass[8pt]{article}\n"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\usepackage[a4paper, total={8in, 10in}]{geometry}"
    latexfile+="\\begin{document}\n"
    latexfile+= "\\centering\n"


    for k in bm.keys():
        latexfile+=f"{k}th Boundary Map's kernel is spanned by\n "
        latexfile+="\\begin{align*}"
        count=0

        for b in bm[k].right_kernel().basis_matrix().rows():

            latexfile+=latex(b.column()) +",\n"
            count+=1
            if count==10:
                latexfile += "\\end{align*}\n"
                latexfile+="\\newline"
                latexfile+="\\begin{align*}\n"

                count=0
                
        latexfile += "\\end{align*}\n"
    latexfile+="\\end{document}"
    file_name=f"So{n}BndryKernel.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 
def printBoundaryMap(n):
    bm=getBndryMaps(n)
    latexfile="\\documentclass[8pt]{article}\n"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\usepackage[a4paper, total={8in, 10in}]{geometry}"
    latexfile+="\\begin{document}\n"
    latexfile+= "\\centering\n"
    for k in bm.keys():

        latexfile+=f"{n*(n-1)//2 + 1 - k}th Boundary Map is given by\n "
        latexfile+="\\begin{align*}"

        latexfile+= latex(bm[k])
        latexfile += "\\end{align*}\n"

    latexfile+="\\end{document}"
    file_name=f"So{n}BndryMap.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 
def printBndryBySignBlock(n,Basis=None):

    bm=getBndryMapByBlocks(n,Basis=Basis) 
    latexfile="\\documentclass[8pt]{article}\n"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\usepackage[a4paper, total={8in, 10in}]{geometry}"
    latexfile+="\\begin{document}\n"
    latexfile+= "\\centering\n"
    for k in bm.keys():

        latexfile+=f"{k}th Boundary Map is given by\n "
        latexfile+="\\begin{align*}"

        latexfile+= latex(bm[k])
        latexfile += "\\end{align*}\n"

    latexfile+="\\end{document}"
    file_name=f"So{n}BndryMapByBlockWithChangeOfbasis.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 
Basism=matrix([[1,0,0,0],[1,1,0,0],[0,1,1,0],[0,0,1,1]],base_ring=GF(Integer(2)))
def printSoChainAndHome(n,showmatrix=False,Basis=None):
    latexfile="\\documentclass[8pt]{article}"
    latexfile+="\\usepackage{amsmath}\n"
    latexfile+="\\usepackage{amssymb}\n"
    latexfile+="\\newcommand{\\Z}{\\mathbb{Z}}\n"

    latexfile+="\\begin{document}\n"

    latexfile+=(f" \\( So( {n} )\\)\n")
    
    bmsage=getBndryMapByBlocks(n,Basis=Basis) 


    C=ChainComplex(bmsage,base_ring=GF(Integer(2)))
    Hom= C.homology()
    Betti=C.betti()
    latexfile+="\\begin{align*}\n"

    latexfile+=latex(C)
    latexfile += "\\end{align*}\n"
    latexfile+="\\begin{align*}\n"

    latexfile+=latex(Hom)
    latexfile += "\\end{align*}\n"
    
    latexfile+="\\end{document}"
    file_name = f"So{n}.tex"
    with open(file_name,"w") as file:
        file.write(latexfile)
    subprocess.run(["pdflatex", file_name]) 
#print(Basism)
#printSignTransBlock(3,Basis=Basism)       
            

#for n in range(3,5):
#    printImageOfBndry(n)    
#    printKernelOfBndry(n)    
#printSo(3)
printSoChainAndHome(3)

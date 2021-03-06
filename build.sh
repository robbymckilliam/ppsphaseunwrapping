#batch file to compile figures and latex
cd plots 
mpost -interaction=nonstopmode circstatfigcube.mp 
mpost -interaction=nonstopmode circstatfigquad.mp 
mpost -interaction=nonstopmode circstatfig.mp 
mpost -interaction=nonstopmode circstatfigzero.mp 
mpost -interaction=nonstopmode tesselationfigures.mp 
cd .. 

cd code/data 
mpost -interaction=nonstopmode plot3.mp 
mpost -interaction=nonstopmode plot5.mp 
mpost -interaction=nonstopmode bench.mp 
cd ../../ 

pdflatex -shell-escape paper.tex 
bibtex paper.aux 
pdflatex -shell-escape paper.tex 
pdflatex -shell-escape paper.tex 

#pdflatex paper2.tex 
#bibtex paper2.aux 
#pdflatex paper2.tex 
#pdflatex paper2.tex 
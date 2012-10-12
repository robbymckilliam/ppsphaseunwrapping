#batch file to compile figures and latex

cd plots 
mpost -interaction=nonstopmode circstatfigcube.mp 
mpost -interaction=nonstopmode circstatfigquad.mp 
mpost -interaction=nonstopmode circstatfig.mp 
mpost -interaction=nonstopmode circstatfigzero.mp 
mpost -interaction=nonstopmode tesselationfigures.mp 
mpost -interaction=nonstopmode 
distplots.class.distributions.circular.ProjectedNormalDistribution.var3.0.mp 
cd .. 

cd code 
mpost -interaction=nonstopmode gaussianplot4.mp 
mpost -interaction=nonstopmode lseplot.mp 
cd .. 

pdflatex paper.tex 
bibtex paper.aux 
pdflatex paper.tex 
pdflatex paper.tex 

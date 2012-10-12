#batch file to compile all the figures from 
#my thesis. If you haven't go rubypost
#installed you will want to uncomment
#the gem line below

#gem install rubypost

cd plots 
ruby1.9.1 circstatfigcube.rb
ruby1.9.1 circstatfigquad.rb
ruby1.9.1 circstatfig.rb
ruby1.9.1 circstatfigzero.rb
ruby1.9.1 tesselationfigures.rb
mpost -interaction=nonstopmode distplots.class.distributions.circular.ProjectedNormalDistribution.var3.0.mp
cd .. 

cd code 
mpost -interaction=nonstopmode gaussianplot4.mp 
mpost -interaction=nonstopmode lseplot.mp 
cd .. 

pdflatex paper.tex
bibtex paper.aux
pdflatex paper.tex
pdflatex paper.tex

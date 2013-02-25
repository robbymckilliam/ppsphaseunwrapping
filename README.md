Repository containing code for the paper 

R. McKilliam, B. Quinn, I. V. L. Clarkson, B. Moran and B. Vellambi. '"Polynomial phase estimation by phase unwrapping"

The main latex file is paper.tex.

The bash script 

build.sh 

with build the figures and compile the latex file to a pdf using pdflatex. The code/ folder contains the output data and software for generating the plots in the paper.  The bash script

runsim.sh 

will run the Monte-Carlo simualtions.  The bash script

benchmark.sh

will run the benchmarks.  You need a working Java Virtual Machine and scala2.9.x or some compatible version to run this.  This software makes use of nearest lattices point algorithms contained in the java library available at 

https://github.com/robbymckilliam/pubsim.  

An executable for this library is packaged inside this repository.

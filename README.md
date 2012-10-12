Repository containing code for the paper 

R. McKilliam, B. Quinn, I. V. L. Clarkson, B. Moran and B. Vellambi. '"Polynomial phase estimation by phase unwrapping"

The main latex file is paper.tex.

The script build.sh with build the figures and compile the latex file to a pdf using pdflatex.  I noticed the script having problem calling bibtex, so you may need to run bibtex paper.aux yourself.

The code folder contains the output data and scripts for generating the plots in the paper.  An existing java library was used for the Monte Carlo simlations.  Execute the script 

runsim.sh 

You need a working jvm and scala2.9.x or some compatible version to run this.  These scripts make use of nearest point algorithms for An* (and other lattices) contained in the java library available at 

https://github.com/robbymckilliam/pubsim.  

An executable for this library is packaged inside this repository.

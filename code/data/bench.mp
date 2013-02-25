prologues := 2;
filenametemplate "%j-%c.mps";

input graph;
input sarith;

%input TEX; 
%TEXPRE("%&latex" & char(10) & 
%"\documentclass[12pt]{article}\begin{document}"); 
%TEXPOST("\end{document}"); 

verbatimtex
%&latex
\documentclass[8pt]{extarticle}
  \begin{document}
  etex

init_numbers(btex $-$ etex,btex $1$ etex, btex ${\times}10$ etex, btex ${}^-$ etex, btex ${}^2$ etex);

picture ORIGINAL_PICTURE;
ORIGINAL_PICTURE := currentpicture;

picture dd; dd :=
dashpattern(on 4bp off 2bp on 1bp off 2bp);

picture smalldash; smalldash :=
dashpattern(on 2bp off 2bp);

%plot pictures
picture diamond;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm rotated 45;
diamond := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture square;
draw unitsquare shifted(-0.5,-0.5) scaled 0.14cm rotated 45;
square := currentpicture; currentpicture := ORIGINAL_PICTURE;
  
picture circle; 
draw fullcircle scaled 0.17cm;
circle := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture fillcirc;
fill fullcircle scaled 0.11cm;
fillcirc := currentpicture; currentpicture := ORIGINAL_PICTURE;
 
picture cross; 
draw (-0.09cm,0.0)--(0.09cm,0.0);
draw (0,-0.09cm)--(0,0.09cm);
cross := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture times; times := cross rotated 45;
%%%%

%legend
picture legendone; 
sx = 0.2cm;
sh = 0.2cm;
ups = 0.35cm;

trans := ups*0;
label.rt(btex LSU etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) dashed evenly;

trans := ups*1;
label.rt(btex LSU (theory) etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans);

legendone := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 6cm;
plotwidth = 7.6cm;


beginfig(1);

draw begingraph(plotwidth, plotheight) ;
  setcoords(log,log);
  glabel.bot(btex SNR (dB) etex, OUT) shifted(0cm, 1cm);
  glabel.lft(btex MSE $\mu_3$ etex, OUT) rotated 90 shifted(1.5cm, 1.8cm);

  gdraw "benchMbest3snr0";
  gdraw "benchMbest3snr10";
  gdraw "benchMbest3snr20";
  gdraw "benchHAF3" dashed evenly;
    
 draw legendone shifted(4.6cm, plotheight-2.2cm);
 
endgraph;

endfig;

end;


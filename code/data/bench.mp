prologues := 2;
filenametemplate "%j-%c.mps";

input graph;
input sarith;

%input TEX; 
%TEXPRE("%&latex" & char(10) & 
%"\documentclass[12pt]{article}\usepackage{units}\begin{document}");
%TEXPOST("\end{document}"); 

verbatimtex
%&latex
\documentclass[8pt]{extarticle}
\usepackage{units}
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
label.rt(btex HAF etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) dashed evenly;

trans := ups*1;
label.rt(btex LSU etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans);

legendone := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 6cm;
plotwidth = 7.6cm;


beginfig(1);

draw begingraph(plotwidth, plotheight) ;
  setcoords(log,log);
  glabel.bot(btex $N$ etex, OUT) shifted(1.05cm, 1cm);
  glabel.lft(btex seconds etex, OUT) rotated 90 shifted(1.3cm, 1.0cm);

  gdraw "benchMbest3snr0"; glabel.ulft(btex \unit[0]{dB} etex, 10);
  gdraw "benchMbest3snr5"; glabel.lrt(btex \unit[5]{dB} etex, 7);
  gdraw "benchMbest3snr10"; glabel.ulft(btex \unit[10]{dB} etex, 7);
  gdraw "benchMbest3snr20"; glabel.lrt(btex \unit[20]{dB} etex, 10);
  gdraw "benchHAF3" dashed evenly;

  %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  %otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $100$ etex, 100);
  otick.lft(btex $1$ etex, "1");
  otick.lft(btex $10^{-2}$ etex, "1e-2");
  otick.lft(btex $10^{-4}$ etex, "1e-4");
  autogrid(otick.bot,);
  
 draw legendone shifted(0.1cm, plotheight-0.85cm);
 
endgraph;

endfig;

verbatimtex
\end{document}
etex

end;


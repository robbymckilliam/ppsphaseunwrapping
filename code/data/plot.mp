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

trans := ups*4;
label.rt(btex CRB etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans);

trans := ups*3;
label.rt(btex K-Best etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw fillcirc shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*6;
label.rt(btex Babai etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw circle shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

%trans := 0;
%label.rt(btex Sphere decoder etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
%draw circle shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*1;
label.rt(btex LSU (theory) etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) dashed evenly;

trans := ups*5;
label.rt(btex DPT etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw cross shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*2;
label.rt(btex CPF etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw square shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

legendone := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 6cm;
plotwidth = 7.6cm;


beginfig(1);

draw begingraph(plotwidth, plotheight) ;
  setcoords(linear,log);
  glabel.bot(btex SNR etex, OUT) shifted(0cm, 1cm);
  glabel.lft(btex MSE etex, OUT) rotated 90 shifted(1.5cm, 1.6cm);

  %gdraw "BabaiEstimatorN10m3Gaussianp3" plot(times);
  %gdraw "DPTEstimatorN39m3p3" plot(cross);
  %gdraw "CubicPhaseFunctionN39m3p3" plot(square);
  %gdraw "MbestEstimatorN39m3p3" plot(fillcirc);
  %gdraw "SphereDecoderEstimatorN39m3p3" plot(circle);

  %gdraw "GaussCRBN39m3p3"; glabel.lrt(btex $N=39$ etex, 0);
  %gdraw "LSUasympGaussianN39m3p3" dashed evenly; 

  gdraw "BabaiEstimatorN199m3p3" plot(circle);
  gdraw "DPTEstimatorN199m3p3" plot(cross);
  gdraw "CubicPhaseFunctionN199m3p3" plot(square);
  gdraw "MbestEstimatorN199m3p3" plot(fillcirc);
  %gdraw "SphereDecoderEstimatorN199m3p3" plot(circle);

  gdraw "GaussCRBN199m3p3"; %glabel.lrt(btex $N=199$ etex, 0);
  gdraw "LSUasympGaussianN199m3p3" dashed evenly; 
  
  %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-4}$ etex, 0.0001);
  otick.lft(btex $10^{-6}$ etex, "1e-6");
  otick.lft(btex $10^{-8}$ etex, "1e-8");
  otick.lft(btex $10^{-10}$ etex, "1e-10");
  otick.lft(btex $10^{-12}$ etex, "1e-12");
  otick.lft(btex $10^{-14}$ etex, "1e-14");
  %otick.lft(btex $10^{-16}$ etex, "1e-16");
  autogrid(otick.bot,);
  
 draw legendone shifted(4.6cm, plotheight-2.2cm);
 
endgraph

endfig;

beginfig(2);

draw begingraph(plotwidth, plotheight) ;
  setcoords(linear,log);
  glabel.bot(btex SNR etex, OUT) shifted(0cm, 1cm);
  glabel.lft(btex MSE etex, OUT) rotated 90 shifted(1.5cm, 1.6cm);

  %gdraw "BabaiEstimatorN10m3Gaussianp3" plot(times);
  %gdraw "DPTEstimatorN39m3p3" plot(cross);
  %gdraw "CubicPhaseFunctionN39m3p3" plot(square);
  %gdraw "MbestEstimatorN39m3p3" plot(fillcirc);
  %gdraw "SphereDecoderEstimatorN39m3p3" plot(circle);

  %gdraw "GaussCRBN39m3p3"; glabel.lrt(btex $N=39$ etex, 0);
  %gdraw "LSUasympGaussianN39m3p3" dashed evenly; 

  gdraw "BabaiEstimatorN199m4p4" plot(circle);
  gdraw "DPTEstimatorN199m4p4" plot(cross);
  %gdraw "CubicPhaseFunctionN199m3p3" plot(square);
  gdraw "MbestEstimatorN199m4p4" plot(fillcirc);
  %gdraw "SphereDecoderEstimatorN199m3p3" plot(circle);

  gdraw "GaussCRBN199m4p4"; %glabel.lrt(btex $N=199$ etex, 0);
  gdraw "LSUasympGaussianN199m4p4" dashed evenly; 
  
  setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  %otick.lft(btex $0.01$ etex, 0.01);
  %otick.lft(btex $10^{-4}$ etex, 0.0001);
  %otick.lft(btex $10^{-6}$ etex, "1e-6");
  %otick.lft(btex $10^{-8}$ etex, "1e-8");
  %otick.lft(btex $10^{-10}$ etex, "1e-10");
  %otick.lft(btex $10^{-12}$ etex, "1e-12");
  %otick.lft(btex $10^{-14}$ etex, "1e-14");
  %otick.lft(btex $10^{-16}$ etex, "1e-16");
  %autogrid(otick.bot,);
  
 draw legendone shifted(4.6cm, plotheight-2.2cm);
 
endgraph

endfig;



end;


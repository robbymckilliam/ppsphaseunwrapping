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
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm;
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

trans := ups*1;
label.rt(btex LSE etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw diamond shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := 0;
label.rt(btex CRB etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans);

trans := ups*3;
label.rt(btex K-Best etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw fillcirc shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*2;
label.rt(btex Babai etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw times shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*5;
label.rt(btex Sphere decoder etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw circle shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*4;
label.rt(btex LSU (theory) etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) dashed evenly;

legendone := currentpicture; currentpicture := ORIGINAL_PICTURE;

%plotheight = 9.7cm;
%plotwidth = 7.6cm;
plotheight = 7cm;
plotwidth = 5cm;

%%%%%%%
% m = 2 Gaussian errors plot
beginfig(1);

draw begingraph(plotwidth, plotheight) ;
  setcoords(log,log);
  glabel.bot(btex $\sigma_c^2$ etex, OUT) shifted(1cm, 1.2cm);
  glabel.lft(btex MSE $\mu_2$ etex, OUT) rotated 90 shifted(1.5cm, 0.2cm);
 
  gdraw "BabaiEstimatorN10m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "DPTEstimatorN10m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "KitchenEstimatorN10m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "MaximumLikelihoodN10m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "MbestEstimatorN10m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "SphereDecoderEstimatorN10m2Gaussianp2" withcolor(0.6,0.6,0.6);

  gdraw "BabaiEstimatorN50m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "DPTEstimatorN50m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "KitchenEstimatorN50m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "MaximumLikelihoodN50m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "MbestEstimatorN50m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "SphereDecoderEstimatorN50m2Gaussianp2" withcolor(0.6,0.6,0.6);

  gdraw "BabaiEstimatorN200m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "DPTEstimatorN200m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "KitchenEstimatorN200m2Gaussianp2" withcolor(0.6,0.6,0.6);
%gdraw "MaximumLikelihoodN200m2Gaussianp2" withcolor(0.6,0.6,0.6);
  gdraw "MbestEstimatorN200m2Gaussianp2" withcolor(0.6,0.6,0.6);
  %gdraw "SphereDecoderEstimatorN200m2Gaussianp2" withcolor(0.6,0.6,0.6);

  gdraw "BabaiEstimatorN10m2Gaussianp2" plot(times);
%gdraw "DPTEstimatorN10m2Gaussianp2" 
  %gdraw "KitchenEstimatorN10m2Gaussianp2" plot(square);
  gdraw "MaximumLikelihoodN10m2Gaussianp2" plot(diamond);
  gdraw "MbestEstimatorN10m2Gaussianp2" plot(fillcirc);
  gdraw "SphereDecoderEstimatorN10m2Gaussianp2" plot(circle);

  gdraw "BabaiEstimatorN50m2Gaussianp2" plot(times);
%gdraw "DPTEstimatorN50m2Gaussianp2" 
  %gdraw "KitchenEstimatorN50m2Gaussianp2" plot(square);
  gdraw "MaximumLikelihoodN50m2Gaussianp2" plot(diamond);
  gdraw "MbestEstimatorN50m2Gaussianp2" plot(fillcirc);
  gdraw "SphereDecoderEstimatorN50m2Gaussianp2" plot(circle);

  gdraw "BabaiEstimatorN200m2Gaussianp2" plot(times);
%gdraw "DPTEstimatorN50m2Gaussianp2" 
  %gdraw "KitchenEstimatorN200m2Gaussianp2" plot(square);
  %gdraw "MaximumLikelihoodN200m2Gaussianp2" plot(diamond);
  gdraw "MbestEstimatorN200m2Gaussianp2" plot(fillcirc);
  %gdraw "SphereDecoderEstimatorN200m2Gaussianp2" plot(circle);

  gdraw "GaussCRBN10m2p2"; glabel.lrt(btex $N=10$ etex, 0) ;
  gdraw "GaussCRBN50m2p2"; glabel.lrt(btex $N=50$ etex, 0) ;
  gdraw "GaussCRBN200m2p2"; glabel.lrt(btex $N=200$ etex, 0) ;

  gdraw "LSUasympGaussianN10m2p2" dashed evenly; 
  gdraw "LSUasympGaussianN50m2p2" dashed evenly;
  gdraw "LSUasympGaussianN200m2p2" dashed evenly;

   %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-4}$ etex, 0.0001);
  otick.lft(btex $10^{-6}$ etex, "1e-6");
  otick.lft(btex $10^{-8}$ etex, "1e-8");
  otick.lft(btex $10^{-10}$ etex, "1e-10");
  otick.lft(btex $10^{-12}$ etex, "1e-12");
  autogrid(otick.bot,);
  
 draw legendone shifted(0.05cm, plotheight - 2.2cm);
 
endgraph;

draw begingraph(plotwidth, plotheight) ;
  setcoords(log,log);
  glabel.bot(btex $\sigma^2$ etex, OUT) shifted(1cm, 1.2cm);
  glabel.lft(btex MSE $\mu_2$ etex, OUT) rotated 90 shifted(1.5cm, 0.5cm);
  
  gdraw "BabaiEstimatorN10m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "DPTEstimatorN10m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "KitchenEstimatorN10m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "MaximumLikelihoodN10m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "MbestEstimatorN10m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "SphereDecoderEstimatorN10m2WUp2" withcolor(0.6,0.6,0.6);

  gdraw "BabaiEstimatorN50m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "DPTEstimatorN50m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "KitchenEstimatorN50m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "MaximumLikelihoodN50m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "MbestEstimatorN50m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "SphereDecoderEstimatorN50m2WUp2" withcolor(0.6,0.6,0.6);

  gdraw "BabaiEstimatorN200m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "DPTEstimatorN200m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "KitchenEstimatorN200m2WUp2" withcolor(0.6,0.6,0.6);
%gdraw "MaximumLikelihoodN200m2WUp2" withcolor(0.6,0.6,0.6);
  gdraw "MbestEstimatorN200m2WUp2" withcolor(0.6,0.6,0.6);
  %gdraw "SphereDecoderEstimatorN200m2WUp2" withcolor(0.6,0.6,0.6);

  gdraw "BabaiEstimatorN10m2WUp2" plot(times);
%gdraw "DPTEstimatorN10m2WUp2" 
  %gdraw "KitchenEstimatorN10m2WUp2" plot(square);
  gdraw "MaximumLikelihoodN10m2WUp2" plot(diamond);
  gdraw "MbestEstimatorN10m2WUp2" plot(fillcirc);
  gdraw "SphereDecoderEstimatorN10m2WUp2" plot(circle);

  gdraw "BabaiEstimatorN50m2WUp2" plot(times);
%gdraw "DPTEstimatorN50m2WUp2" 
  %gdraw "KitchenEstimatorN50m2WUp2" plot(square);
  gdraw "MaximumLikelihoodN50m2WUp2" plot(diamond);
  gdraw "MbestEstimatorN50m2WUp2" plot(fillcirc);
  gdraw "SphereDecoderEstimatorN50m2WUp2" plot(circle);

  gdraw "BabaiEstimatorN200m2WUp2" plot(times);
%gdraw "DPTEstimatorN50m2WUp2" 
  %gdraw "KitchenEstimatorN200m2WUp2" plot(square);
  %gdraw "MaximumLikelihoodN200m2WUp2" plot(diamond);
  gdraw "MbestEstimatorN200m2WUp2" plot(fillcirc);
  %gdraw "SphereDecoderEstimatorN200m2WUp2" plot(circle);

  gdraw "LSUasympWUN10m2p2" dashed evenly; glabel.lrt(btex $N=10$ etex, 19); 
  gdraw "LSUasympWUN50m2p2" dashed evenly; glabel.lrt(btex $N=50$ etex, 19);
  gdraw "LSUasympWUN200m2p2" dashed evenly; glabel.lrt(btex $N=200$ etex, 19);

   %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-4}$ etex, 0.0001);
  otick.lft(btex $10^{-6}$ etex, "1e-6");
  otick.lft(btex $10^{-8}$ etex, "1e-8");
  otick.lft(btex $10^{-10}$ etex, "1e-10");
  otick.lft(btex $10^{-12}$ etex, "1e-12");
  autogrid(otick.bot,);
  
 %draw legendone shifted(0.05cm, 7.2cm);
 
endgraph shifted(plotwidth + 1.5cm, 0);

endfig;

end;

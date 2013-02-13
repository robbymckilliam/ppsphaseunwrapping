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

trans := ups*2;
label.rt(btex CRB etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans);

trans := ups*3;
label.rt(btex LSU etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw fillcirc shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

% trans := ups*6;
% label.rt(btex Babai etex, (1cm,sh + trans));
% %draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
% draw circle shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

%trans := 0;
%label.rt(btex Sphere decoder etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
%draw circle shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*1;
label.rt(btex LSU (theory) etex, (1cm,sh + trans));
draw (sx,sh + trans)--(1cm,sh + trans) dashed evenly;

trans := ups*5;
label.rt(btex HAF etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw cross shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

trans := ups*4;
label.rt(btex CPF etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw square shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

legendone := currentpicture; currentpicture := ORIGINAL_PICTURE;


%legend
picture legendtwo; 
sx = 0.2cm;
sh = 0.2cm;
ups = 0.35cm;

%trans := ups*2;
%label.rt(btex CRB etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans);

trans := ups*3;
label.rt(btex LSU etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw fillcirc shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

%trans := ups*1;
%label.rt(btex LSU (theory) etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) dashed evenly;

trans := ups*4;
label.rt(btex Zhou Wang etex, (1cm,sh + trans));
%draw (sx,sh + trans)--(1cm,sh + trans) withcolor(0.6,0.6,0.6);
draw cross shifted 1/2[(sx,sh + trans),(1cm,sh + trans)];

legendtwo := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 6cm;
plotwidth = 7.6cm;


beginfig(1);

draw begingraph(plotwidth, plotheight) ;
  setcoords(linear,log);
  glabel.bot(btex SNR (dB) etex, OUT) shifted(0cm, 1cm);
  glabel.lft(btex MSE $\mu_3$ etex, OUT) rotated 90 shifted(1.5cm, 1.8cm);

  %gdraw "BabaiN10m3Gaussianp3" plot(times);
  %gdraw "HAFN39m3p3" plot(cross);
  %gdraw "CPFN39m3p3" plot(square);
  %gdraw "MbestN39m3p3" plot(fillcirc);
  %gdraw "SphereDecoderN39m3p3" plot(circle);

  %gdraw "GaussCRBN39m3p3"; glabel.lrt(btex $N=39$ etex, 0);
  %gdraw "LSUasympGaussianN39m3p3" dashed evenly; 

  %gdraw "Babaim3smallp3" plot(circle);
  gdraw "HAFm3smallp3" plot(cross);
  gdraw "CPFm3smallp3" plot(square);
  gdraw "Mbestm3smallp3" plot(fillcirc);
  %gdraw "ZWN199m3p3small";
  %gdraw "SphereDecoderN199m3p3" plot(circle);

  gdraw "crbm3p3"; %glabel.lrt(btex $N=199$ etex, 0);
  gdraw "lsucltm3p3" dashed evenly; 
  
  %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  %otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-3}$ etex, 0.001);
  otick.lft(btex $10^{-5}$ etex, "1e-5");
  otick.lft(btex $10^{-7}$ etex, "1e-7");
  otick.lft(btex $10^{-9}$ etex, "1e-9");
  otick.lft(btex $10^{-11}$ etex, "1e-11");
  otick.lft(btex $10^{-13}$ etex, "1e-13");
  otick.lft(btex $10^{-15}$ etex, "1e-15");
  autogrid(otick.bot,);
  
 draw legendone shifted(4.6cm, plotheight-2.2cm);
 
endgraph;

endfig;


beginfig(2);

draw begingraph(plotwidth, plotheight) ;
  setcoords(linear,log);
  glabel.bot(btex SNR (dB) etex, OUT) shifted(0cm, 1cm);
  glabel.lft(btex MSE $\mu_3$ etex, OUT) rotated 90 shifted(1.5cm, 1.6cm);

  %gdraw "BabaiN10m3Gaussianp3" plot(times);
  %gdraw "HAFN39m3p3" plot(cross);
  %gdraw "CPFN39m3p3" plot(square);
  %gdraw "MbestN39m3p3" plot(fillcirc);
  %gdraw "SphereDecoderN39m3p3" plot(circle);

  %gdraw "GaussCRBN39m3p3"; glabel.lrt(btex $N=39$ etex, 0);
  %gdraw "LSUasympGaussianN39m3p3" dashed evenly; 

  %gdraw "Babaim3bigp3" plot(circle);
  gdraw "HAFm3bigp3" plot(cross);
  gdraw "CPFm3bigp3" plot(square);
  gdraw "Mbestm3bigp3" plot(fillcirc);
  %gdraw "ZWN199m3p3big";
  %gdraw "SphereDecoderN199m3p3" plot(circle);

  gdraw "crbm3p3"; %glabel.lrt(btex $N=199$ etex, 0);
  gdraw "lsucltm3p3" dashed evenly; 
  
  %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  %otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-3}$ etex, 0.001);
  otick.lft(btex $10^{-5}$ etex, "1e-5");
  otick.lft(btex $10^{-7}$ etex, "1e-7");
  otick.lft(btex $10^{-9}$ etex, "1e-9");
  otick.lft(btex $10^{-11}$ etex, "1e-11");
  otick.lft(btex $10^{-13}$ etex, "1e-13");
  otick.lft(btex $10^{-15}$ etex, "1e-15");
  autogrid(otick.bot,);
  
 %draw legendone shifted(4.6cm, plotheight-2.7cm);
 
endgraph;

endfig;



beginfig(3);

draw begingraph(plotwidth, plotheight) ;
  setcoords(linear,log);
  glabel.bot(btex SNR (dB)etex, OUT) shifted(0cm, 1cm);
  glabel.lft(btex MSE $\mu_3$ etex, OUT) rotated 90 shifted(1.5cm, 1.6cm);

  %gdraw "BabaiN10m3Gaussianp3" plot(times);
  %gdraw "HAFN39m3p3" plot(cross);
  %gdraw "CPFN39m3p3" plot(square);
  %gdraw "MbestN39m3p3" plot(fillcirc);
  %gdraw "SphereDecoderN39m3p3" plot(circle);

  %gdraw "GaussCRBN39m3p3"; glabel.lrt(btex $N=39$ etex, 0);
  %gdraw "LSUasympGaussianN39m3p3" dashed evenly; 

  %gdraw "Babaim3bigrangep3" plot(circle);
  %gdraw "HAFN199m3p3big" plot(cross);
  %gdraw "CPFN199m3p3big" plot(square);
  gdraw "Mbestm3bigrangep3" plot(fillcirc);
  gdraw "ZWm3bigp3" plot(cross);
  %gdraw "SphereDecoderN199m3p3" plot(circle);

  gdraw "crbm3rangep3"; %glabel.lrt(btex $N=199$ etex, 0);
  gdraw "lsucltm3rangep3" dashed evenly; 
  
  %setrange(whatever,whatever, whatever,whatever);
  %otick.lft(btex $10$ etex, 10);
  %otick.lft(btex $0.01$ etex, 0.01);
  otick.lft(btex $10^{-3}$ etex, 0.001);
  otick.lft(btex $10^{-5}$ etex, "1e-5");
  otick.lft(btex $10^{-7}$ etex, "1e-7");
  otick.lft(btex $10^{-9}$ etex, "1e-9");
  otick.lft(btex $10^{-11}$ etex, "1e-11");
  otick.lft(btex $10^{-13}$ etex, "1e-13");
  otick.lft(btex $10^{-15}$ etex, "1e-15");
  otick.lft(btex $10^{-17}$ etex, "1e-17");
  otick.bot(btex $-5$ etex, -5);
  otick.bot(btex $0$ etex, 0);
  otick.bot(btex $5$ etex, 5);
  otick.bot(btex $10$ etex, 10);
  otick.bot(btex $15$ etex, 15);
  otick.bot(btex $20$ etex, 20);
  otick.bot(btex $20$ etex, 20);
  otick.bot(btex $25$ etex, 25);
  otick.bot(btex $30$ etex, 30);
  otick.bot(btex $35$ etex, 35);
  %autogrid(otick.bot,);
  
draw legendtwo shifted(4.6cm, plotheight-1.9cm);
 
endgraph;

endfig;


end;


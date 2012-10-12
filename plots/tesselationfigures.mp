prologues := 2;
filenametemplate "%j-%c.mps";
verbatimtex
%&latex
\documentclass{minimal}
\begin{document}
etex
picture ORIGINAL_PICTURE;
ORIGINAL_PICTURE := currentpicture;
picture axes;
axes := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture picL;
draw axes ;

fill fullcircle  scaled 0.07cm shifted (-7.2cm, -7.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-6.96cm, -6.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-6.72cm, -4.8cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-6.48cm, -3.5999999999999996cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-6.24cm, -2.4cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-6.0cm, -1.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.76cm, 0.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.52cm, 1.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.28cm, 2.3999999999999995cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.04cm, 3.5999999999999996cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.8cm, 4.8cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-6.0cm, -6.96cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.76cm, -5.76cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.52cm, -4.56cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.279999999999999cm, -3.36cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-5.04cm, -2.16cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.8cm, -0.96cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.56cm, 0.24cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.32cm, 1.44cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.08cm, 2.6399999999999997cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.84cm, 3.84cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.5999999999999996cm, 5.04cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.8cm, -6.72cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.56cm, -5.52cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.319999999999999cm, -4.319999999999999cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-4.08cm, -3.12cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.84cm, -1.92cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.5999999999999996cm, -0.72cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.3599999999999994cm, 0.48cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.1199999999999997cm, 1.68cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.88cm, 2.88cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.6399999999999997cm, 4.08cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.3999999999999995cm, 5.28cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.5999999999999996cm, -6.48cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.36cm, -5.279999999999999cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-3.12cm, -4.08cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.88cm, -2.88cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.6399999999999997cm, -1.68cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.4cm, -0.48cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.16cm, 0.72cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.92cm, 1.92cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.68cm, 3.1199999999999997cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.44cm, 4.32cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.2cm, 5.52cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.4cm, -6.24cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-2.16cm, -5.04cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.92cm, -3.84cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.68cm, -2.6399999999999997cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.44cm, -1.44cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.2cm, -0.24cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.96cm, 0.96cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.72cm, 2.16cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.48cm, 3.3599999999999994cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.24cm, 4.56cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.0cm, 5.76cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-1.2cm, -6.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.96cm, -4.8cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.72cm, -3.5999999999999996cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.48cm, -2.4cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (-0.24cm, -1.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.0cm, 0.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.24cm, 1.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.48cm, 2.4cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.72cm, 3.5999999999999996cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.96cm, 4.8cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.2cm, 6.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.0cm, -5.76cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.24cm, -4.56cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.48cm, -3.3599999999999994cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.72cm, -2.16cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (0.96cm, -0.96cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.2cm, 0.24cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.44cm, 1.44cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.68cm, 2.6399999999999997cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.92cm, 3.84cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.16cm, 5.04cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.4cm, 6.24cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.2cm, -5.52cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.44cm, -4.32cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.68cm, -3.1199999999999997cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (1.92cm, -1.92cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.16cm, -0.72cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.4cm, 0.48cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.6399999999999997cm, 1.68cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.88cm, 2.88cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.12cm, 4.08cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.36cm, 5.279999999999999cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.5999999999999996cm, 6.48cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.3999999999999995cm, -5.28cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.6399999999999997cm, -4.08cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (2.88cm, -2.88cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.1199999999999997cm, -1.68cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.3599999999999994cm, -0.48cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.5999999999999996cm, 0.72cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.84cm, 1.92cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.08cm, 3.12cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.319999999999999cm, 4.319999999999999cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.56cm, 5.52cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.8cm, 6.72cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.5999999999999996cm, -5.04cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (3.84cm, -3.84cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.08cm, -2.6399999999999997cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.32cm, -1.44cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.56cm, -0.24cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.8cm, 0.96cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.04cm, 2.16cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.279999999999999cm, 3.36cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.52cm, 4.56cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.76cm, 5.76cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (6.0cm, 6.96cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (4.8cm, -4.8cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.04cm, -3.5999999999999996cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.28cm, -2.3999999999999995cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.52cm, -1.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (5.76cm, 0.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (6.0cm, 1.2cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (6.24cm, 2.4cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (6.48cm, 3.5999999999999996cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (6.72cm, 4.8cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (6.96cm, 6.0cm) withcolor(0,0,0);

fill fullcircle  scaled 0.07cm shifted (7.2cm, 7.2cm) withcolor(0,0,0);

picL := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture trishaded;
fill ((0, 0)--(1.2cm, 0.24cm)--(0.24cm, 1.2cm)--cycle)  shifted (0.48cm, 0.48cm) withcolor(0.7,0.7,0.7);

fill ((-1.44cm, -1.44cm)--(-2.6399999999999997cm, -1.68cm)--(-1.68cm, -2.6399999999999997cm)--cycle)  shifted (0.48cm, 0.48cm) withcolor(0.7,0.7,0.7);

trishaded := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture trioutlined;
draw ((0, 0)--(1.2cm, 0.24cm)--(0.24cm, 1.2cm)--cycle)  shifted (0.48cm, 0.48cm) withcolor(0.6,0.6,0.6);

draw ((-1.44cm, -1.44cm)--(-2.6399999999999997cm, -1.68cm)--(-1.68cm, -2.6399999999999997cm)--cycle)  shifted (0.48cm, 0.48cm) withcolor(0.6,0.6,0.6);

trioutlined := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture tritesselate;
draw trioutlined  shifted (-7.2cm, -7.2cm);

draw trioutlined  shifted (-6.96cm, -6.0cm);

draw trioutlined  shifted (-6.72cm, -4.8cm);

draw trioutlined  shifted (-6.48cm, -3.5999999999999996cm);

draw trioutlined  shifted (-6.24cm, -2.4cm);

draw trioutlined  shifted (-6.0cm, -1.2cm);

draw trioutlined  shifted (-5.76cm, 0.0cm);

draw trioutlined  shifted (-5.52cm, 1.2cm);

draw trioutlined  shifted (-5.28cm, 2.3999999999999995cm);

draw trioutlined  shifted (-5.04cm, 3.5999999999999996cm);

draw trioutlined  shifted (-4.8cm, 4.8cm);

draw trioutlined  shifted (-6.0cm, -6.96cm);

draw trioutlined  shifted (-5.76cm, -5.76cm);

draw trioutlined  shifted (-5.52cm, -4.56cm);

draw trioutlined  shifted (-5.279999999999999cm, -3.36cm);

draw trioutlined  shifted (-5.04cm, -2.16cm);

draw trioutlined  shifted (-4.8cm, -0.96cm);

draw trioutlined  shifted (-4.56cm, 0.24cm);

draw trioutlined  shifted (-4.32cm, 1.44cm);

draw trioutlined  shifted (-4.08cm, 2.6399999999999997cm);

draw trioutlined  shifted (-3.84cm, 3.84cm);

draw trioutlined  shifted (-3.5999999999999996cm, 5.04cm);

draw trioutlined  shifted (-4.8cm, -6.72cm);

draw trioutlined  shifted (-4.56cm, -5.52cm);

draw trioutlined  shifted (-4.319999999999999cm, -4.319999999999999cm);

draw trioutlined  shifted (-4.08cm, -3.12cm);

draw trioutlined  shifted (-3.84cm, -1.92cm);

draw trioutlined  shifted (-3.5999999999999996cm, -0.72cm);

draw trioutlined  shifted (-3.3599999999999994cm, 0.48cm);

draw trioutlined  shifted (-3.1199999999999997cm, 1.68cm);

draw trioutlined  shifted (-2.88cm, 2.88cm);

draw trioutlined  shifted (-2.6399999999999997cm, 4.08cm);

draw trioutlined  shifted (-2.3999999999999995cm, 5.28cm);

draw trioutlined  shifted (-3.5999999999999996cm, -6.48cm);

draw trioutlined  shifted (-3.36cm, -5.279999999999999cm);

draw trioutlined  shifted (-3.12cm, -4.08cm);

draw trioutlined  shifted (-2.88cm, -2.88cm);

draw trioutlined  shifted (-2.6399999999999997cm, -1.68cm);

draw trioutlined  shifted (-2.4cm, -0.48cm);

draw trioutlined  shifted (-2.16cm, 0.72cm);

draw trioutlined  shifted (-1.92cm, 1.92cm);

draw trioutlined  shifted (-1.68cm, 3.1199999999999997cm);

draw trioutlined  shifted (-1.44cm, 4.32cm);

draw trioutlined  shifted (-1.2cm, 5.52cm);

draw trioutlined  shifted (-2.4cm, -6.24cm);

draw trioutlined  shifted (-2.16cm, -5.04cm);

draw trioutlined  shifted (-1.92cm, -3.84cm);

draw trioutlined  shifted (-1.68cm, -2.6399999999999997cm);

draw trioutlined  shifted (-1.44cm, -1.44cm);

draw trioutlined  shifted (-1.2cm, -0.24cm);

draw trioutlined  shifted (-0.96cm, 0.96cm);

draw trioutlined  shifted (-0.72cm, 2.16cm);

draw trioutlined  shifted (-0.48cm, 3.3599999999999994cm);

draw trioutlined  shifted (-0.24cm, 4.56cm);

draw trioutlined  shifted (0.0cm, 5.76cm);

draw trioutlined  shifted (-1.2cm, -6.0cm);

draw trioutlined  shifted (-0.96cm, -4.8cm);

draw trioutlined  shifted (-0.72cm, -3.5999999999999996cm);

draw trioutlined  shifted (-0.48cm, -2.4cm);

draw trioutlined  shifted (-0.24cm, -1.2cm);

draw trioutlined  shifted (0.0cm, 0.0cm);

draw trioutlined  shifted (0.24cm, 1.2cm);

draw trioutlined  shifted (0.48cm, 2.4cm);

draw trioutlined  shifted (0.72cm, 3.5999999999999996cm);

draw trioutlined  shifted (0.96cm, 4.8cm);

draw trioutlined  shifted (1.2cm, 6.0cm);

draw trioutlined  shifted (0.0cm, -5.76cm);

draw trioutlined  shifted (0.24cm, -4.56cm);

draw trioutlined  shifted (0.48cm, -3.3599999999999994cm);

draw trioutlined  shifted (0.72cm, -2.16cm);

draw trioutlined  shifted (0.96cm, -0.96cm);

draw trioutlined  shifted (1.2cm, 0.24cm);

draw trioutlined  shifted (1.44cm, 1.44cm);

draw trioutlined  shifted (1.68cm, 2.6399999999999997cm);

draw trioutlined  shifted (1.92cm, 3.84cm);

draw trioutlined  shifted (2.16cm, 5.04cm);

draw trioutlined  shifted (2.4cm, 6.24cm);

draw trioutlined  shifted (1.2cm, -5.52cm);

draw trioutlined  shifted (1.44cm, -4.32cm);

draw trioutlined  shifted (1.68cm, -3.1199999999999997cm);

draw trioutlined  shifted (1.92cm, -1.92cm);

draw trioutlined  shifted (2.16cm, -0.72cm);

draw trioutlined  shifted (2.4cm, 0.48cm);

draw trioutlined  shifted (2.6399999999999997cm, 1.68cm);

draw trioutlined  shifted (2.88cm, 2.88cm);

draw trioutlined  shifted (3.12cm, 4.08cm);

draw trioutlined  shifted (3.36cm, 5.279999999999999cm);

draw trioutlined  shifted (3.5999999999999996cm, 6.48cm);

draw trioutlined  shifted (2.3999999999999995cm, -5.28cm);

draw trioutlined  shifted (2.6399999999999997cm, -4.08cm);

draw trioutlined  shifted (2.88cm, -2.88cm);

draw trioutlined  shifted (3.1199999999999997cm, -1.68cm);

draw trioutlined  shifted (3.3599999999999994cm, -0.48cm);

draw trioutlined  shifted (3.5999999999999996cm, 0.72cm);

draw trioutlined  shifted (3.84cm, 1.92cm);

draw trioutlined  shifted (4.08cm, 3.12cm);

draw trioutlined  shifted (4.319999999999999cm, 4.319999999999999cm);

draw trioutlined  shifted (4.56cm, 5.52cm);

draw trioutlined  shifted (4.8cm, 6.72cm);

draw trioutlined  shifted (3.5999999999999996cm, -5.04cm);

draw trioutlined  shifted (3.84cm, -3.84cm);

draw trioutlined  shifted (4.08cm, -2.6399999999999997cm);

draw trioutlined  shifted (4.32cm, -1.44cm);

draw trioutlined  shifted (4.56cm, -0.24cm);

draw trioutlined  shifted (4.8cm, 0.96cm);

draw trioutlined  shifted (5.04cm, 2.16cm);

draw trioutlined  shifted (5.279999999999999cm, 3.36cm);

draw trioutlined  shifted (5.52cm, 4.56cm);

draw trioutlined  shifted (5.76cm, 5.76cm);

draw trioutlined  shifted (6.0cm, 6.96cm);

draw trioutlined  shifted (4.8cm, -4.8cm);

draw trioutlined  shifted (5.04cm, -3.5999999999999996cm);

draw trioutlined  shifted (5.28cm, -2.3999999999999995cm);

draw trioutlined  shifted (5.52cm, -1.2cm);

draw trioutlined  shifted (5.76cm, 0.0cm);

draw trioutlined  shifted (6.0cm, 1.2cm);

draw trioutlined  shifted (6.24cm, 2.4cm);

draw trioutlined  shifted (6.48cm, 3.5999999999999996cm);

draw trioutlined  shifted (6.72cm, 4.8cm);

draw trioutlined  shifted (6.96cm, 6.0cm);

draw trioutlined  shifted (7.2cm, 7.2cm);

tritesselate := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture rectshaded;
fill ((0.407692307692308cm, 0.561538461538462cm)--(-0.592307692307692cm, 0.361538461538461cm)--(-0.407692307692308cm, -0.561538461538462cm)--(0.592307692307692cm, -0.361538461538461cm)--cycle)  scaled 1.2 withcolor(0.7,0.7,0.7);

rectshaded := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture rectoutlined;
draw ((0.407692307692308cm, 0.561538461538462cm)--(-0.592307692307692cm, 0.361538461538461cm)--(-0.407692307692308cm, -0.561538461538462cm)--(0.592307692307692cm, -0.361538461538461cm)--cycle)  scaled 1.2;

rectoutlined := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture recttesselate;
draw rectoutlined  shifted (-7.2cm, -7.2cm);

draw rectoutlined  shifted (-6.96cm, -6.0cm);

draw rectoutlined  shifted (-6.72cm, -4.8cm);

draw rectoutlined  shifted (-6.48cm, -3.5999999999999996cm);

draw rectoutlined  shifted (-6.24cm, -2.4cm);

draw rectoutlined  shifted (-6.0cm, -1.2cm);

draw rectoutlined  shifted (-5.76cm, 0.0cm);

draw rectoutlined  shifted (-5.52cm, 1.2cm);

draw rectoutlined  shifted (-5.28cm, 2.3999999999999995cm);

draw rectoutlined  shifted (-5.04cm, 3.5999999999999996cm);

draw rectoutlined  shifted (-4.8cm, 4.8cm);

draw rectoutlined  shifted (-6.0cm, -6.96cm);

draw rectoutlined  shifted (-5.76cm, -5.76cm);

draw rectoutlined  shifted (-5.52cm, -4.56cm);

draw rectoutlined  shifted (-5.279999999999999cm, -3.36cm);

draw rectoutlined  shifted (-5.04cm, -2.16cm);

draw rectoutlined  shifted (-4.8cm, -0.96cm);

draw rectoutlined  shifted (-4.56cm, 0.24cm);

draw rectoutlined  shifted (-4.32cm, 1.44cm);

draw rectoutlined  shifted (-4.08cm, 2.6399999999999997cm);

draw rectoutlined  shifted (-3.84cm, 3.84cm);

draw rectoutlined  shifted (-3.5999999999999996cm, 5.04cm);

draw rectoutlined  shifted (-4.8cm, -6.72cm);

draw rectoutlined  shifted (-4.56cm, -5.52cm);

draw rectoutlined  shifted (-4.319999999999999cm, -4.319999999999999cm);

draw rectoutlined  shifted (-4.08cm, -3.12cm);

draw rectoutlined  shifted (-3.84cm, -1.92cm);

draw rectoutlined  shifted (-3.5999999999999996cm, -0.72cm);

draw rectoutlined  shifted (-3.3599999999999994cm, 0.48cm);

draw rectoutlined  shifted (-3.1199999999999997cm, 1.68cm);

draw rectoutlined  shifted (-2.88cm, 2.88cm);

draw rectoutlined  shifted (-2.6399999999999997cm, 4.08cm);

draw rectoutlined  shifted (-2.3999999999999995cm, 5.28cm);

draw rectoutlined  shifted (-3.5999999999999996cm, -6.48cm);

draw rectoutlined  shifted (-3.36cm, -5.279999999999999cm);

draw rectoutlined  shifted (-3.12cm, -4.08cm);

draw rectoutlined  shifted (-2.88cm, -2.88cm);

draw rectoutlined  shifted (-2.6399999999999997cm, -1.68cm);

draw rectoutlined  shifted (-2.4cm, -0.48cm);

draw rectoutlined  shifted (-2.16cm, 0.72cm);

draw rectoutlined  shifted (-1.92cm, 1.92cm);

draw rectoutlined  shifted (-1.68cm, 3.1199999999999997cm);

draw rectoutlined  shifted (-1.44cm, 4.32cm);

draw rectoutlined  shifted (-1.2cm, 5.52cm);

draw rectoutlined  shifted (-2.4cm, -6.24cm);

draw rectoutlined  shifted (-2.16cm, -5.04cm);

draw rectoutlined  shifted (-1.92cm, -3.84cm);

draw rectoutlined  shifted (-1.68cm, -2.6399999999999997cm);

draw rectoutlined  shifted (-1.44cm, -1.44cm);

draw rectoutlined  shifted (-1.2cm, -0.24cm);

draw rectoutlined  shifted (-0.96cm, 0.96cm);

draw rectoutlined  shifted (-0.72cm, 2.16cm);

draw rectoutlined  shifted (-0.48cm, 3.3599999999999994cm);

draw rectoutlined  shifted (-0.24cm, 4.56cm);

draw rectoutlined  shifted (0.0cm, 5.76cm);

draw rectoutlined  shifted (-1.2cm, -6.0cm);

draw rectoutlined  shifted (-0.96cm, -4.8cm);

draw rectoutlined  shifted (-0.72cm, -3.5999999999999996cm);

draw rectoutlined  shifted (-0.48cm, -2.4cm);

draw rectoutlined  shifted (-0.24cm, -1.2cm);

draw rectoutlined  shifted (0.0cm, 0.0cm);

draw rectoutlined  shifted (0.24cm, 1.2cm);

draw rectoutlined  shifted (0.48cm, 2.4cm);

draw rectoutlined  shifted (0.72cm, 3.5999999999999996cm);

draw rectoutlined  shifted (0.96cm, 4.8cm);

draw rectoutlined  shifted (1.2cm, 6.0cm);

draw rectoutlined  shifted (0.0cm, -5.76cm);

draw rectoutlined  shifted (0.24cm, -4.56cm);

draw rectoutlined  shifted (0.48cm, -3.3599999999999994cm);

draw rectoutlined  shifted (0.72cm, -2.16cm);

draw rectoutlined  shifted (0.96cm, -0.96cm);

draw rectoutlined  shifted (1.2cm, 0.24cm);

draw rectoutlined  shifted (1.44cm, 1.44cm);

draw rectoutlined  shifted (1.68cm, 2.6399999999999997cm);

draw rectoutlined  shifted (1.92cm, 3.84cm);

draw rectoutlined  shifted (2.16cm, 5.04cm);

draw rectoutlined  shifted (2.4cm, 6.24cm);

draw rectoutlined  shifted (1.2cm, -5.52cm);

draw rectoutlined  shifted (1.44cm, -4.32cm);

draw rectoutlined  shifted (1.68cm, -3.1199999999999997cm);

draw rectoutlined  shifted (1.92cm, -1.92cm);

draw rectoutlined  shifted (2.16cm, -0.72cm);

draw rectoutlined  shifted (2.4cm, 0.48cm);

draw rectoutlined  shifted (2.6399999999999997cm, 1.68cm);

draw rectoutlined  shifted (2.88cm, 2.88cm);

draw rectoutlined  shifted (3.12cm, 4.08cm);

draw rectoutlined  shifted (3.36cm, 5.279999999999999cm);

draw rectoutlined  shifted (3.5999999999999996cm, 6.48cm);

draw rectoutlined  shifted (2.3999999999999995cm, -5.28cm);

draw rectoutlined  shifted (2.6399999999999997cm, -4.08cm);

draw rectoutlined  shifted (2.88cm, -2.88cm);

draw rectoutlined  shifted (3.1199999999999997cm, -1.68cm);

draw rectoutlined  shifted (3.3599999999999994cm, -0.48cm);

draw rectoutlined  shifted (3.5999999999999996cm, 0.72cm);

draw rectoutlined  shifted (3.84cm, 1.92cm);

draw rectoutlined  shifted (4.08cm, 3.12cm);

draw rectoutlined  shifted (4.319999999999999cm, 4.319999999999999cm);

draw rectoutlined  shifted (4.56cm, 5.52cm);

draw rectoutlined  shifted (4.8cm, 6.72cm);

draw rectoutlined  shifted (3.5999999999999996cm, -5.04cm);

draw rectoutlined  shifted (3.84cm, -3.84cm);

draw rectoutlined  shifted (4.08cm, -2.6399999999999997cm);

draw rectoutlined  shifted (4.32cm, -1.44cm);

draw rectoutlined  shifted (4.56cm, -0.24cm);

draw rectoutlined  shifted (4.8cm, 0.96cm);

draw rectoutlined  shifted (5.04cm, 2.16cm);

draw rectoutlined  shifted (5.279999999999999cm, 3.36cm);

draw rectoutlined  shifted (5.52cm, 4.56cm);

draw rectoutlined  shifted (5.76cm, 5.76cm);

draw rectoutlined  shifted (6.0cm, 6.96cm);

draw rectoutlined  shifted (4.8cm, -4.8cm);

draw rectoutlined  shifted (5.04cm, -3.5999999999999996cm);

draw rectoutlined  shifted (5.28cm, -2.3999999999999995cm);

draw rectoutlined  shifted (5.52cm, -1.2cm);

draw rectoutlined  shifted (5.76cm, 0.0cm);

draw rectoutlined  shifted (6.0cm, 1.2cm);

draw rectoutlined  shifted (6.24cm, 2.4cm);

draw rectoutlined  shifted (6.48cm, 3.5999999999999996cm);

draw rectoutlined  shifted (6.72cm, 4.8cm);

draw rectoutlined  shifted (6.96cm, 6.0cm);

draw rectoutlined  shifted (7.2cm, 7.2cm);

recttesselate := currentpicture; currentpicture := ORIGINAL_PICTURE;

beginfig(1);
clip picL to ((-0.5, -0.5)--(-0.5, 0.5)--(0.5, 0.5)--(0.5, -0.5)--cycle)  scaled 7.93cm;

clip tritesselate to ((-0.5, -0.5)--(-0.5, 0.5)--(0.5, 0.5)--(0.5, -0.5)--cycle)  scaled 7.93cm;

draw trishaded ;

draw picL ;

draw tritesselate  shifted (8.4cm, 0);

draw picL  shifted (8.4cm, 0);

endfig;

beginfig(2);
draw rectshaded ;

clip picL to (-4cm,-3cm)--(4cm,-3cm)--(4cm,3cm)--(-4cm,3cm)--cycle ;

clip recttesselate to (-4cm,-3cm)--(4cm,-3cm)--(4cm,3cm)--(-4cm,3cm)--cycle ;

draw picL ;

draw recttesselate  shifted (0, 0);

draw picL  shifted (0, 0);

endfig;

end;

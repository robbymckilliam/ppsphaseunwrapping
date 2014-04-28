require 'rubypost'

require 'redefine_float_to_s_for_metapost.rb'
require 'matrix'
require 'mathn'

include RubyPost

file = RubyPost::File.new('tesselationfigures')

#create an axis peicture
axes = Picture.new('axes')

#axes.add_drawable('draw (0.15cm, 0)--(-0.15cm, 0);')
#axes.add_drawable('draw (0, 0.15cm)--(0, -0.15cm);')

#~ xaxis = Path.new
#~ xaxis.add_pair(Pair.new(-3.95.cm, 0))
#~ xaxis.add_pair(Pair.new(3.95.cm, 0))
#~ yaxis = Path.new
#~ yaxis.add_pair(Pair.new(0, -3.95.cm))
#~ yaxis.add_pair(Pair.new(0, 3.95.cm))

#~ axes.add_drawable(DoubleArrow.new(xaxis))
#~ axes.add_drawable(DoubleArrow.new(yaxis))
#axes.add_drawable(Label.new(latex("$y$"), Pair.new(0,3.95.cm)).bottom_right)
#axes.add_drawable(Label.new(latex("$x$"), Pair.new(3.95.cm, 0)).top_left)



pic = Picture.new('picL')
pic.add_drawable(Draw.new(axes))

matscale = 1.2

M = Matrix[ [1,0.2], [0.2, 1] ]*(matscale)
(-5..5).each do |x|
	(-5..5).each do |y|
		v = M*Vector[x,y]
		pic.add_drawable(Fill.new(Circle.new()).scale(0.07.cm).translate(v[0].cm,v[1].cm).colour(0,0,0))
	end
end

#####################
# Figure with disconnected triangular tesselating region
#####################
v10 = M*Vector[1,0]
v01 = M*Vector[0,1]
v10n = M*Vector[-2,-1]
v01n = M*Vector[-1,-2]
v11n = M*Vector[-1,-1]
v00 = M*Vector[0,0]
tricyc1 = Path.new
tricyc1.add_pair(Pair.new(0,0))
tricyc1.add_pair(Pair.new(v10[0].cm,v10[1].cm))
tricyc1.add_pair(Pair.new(v01[0].cm,v01[1].cm))
tricyc1.add_pair('cycle')
tricyc2 = Path.new
tricyc2.add_pair(Pair.new(v11n[0].cm,v11n[0].cm))
tricyc2.add_pair(Pair.new(v10n[0].cm,v10n[1].cm))
tricyc2.add_pair(Pair.new(v01n[0].cm,v01n[1].cm))
tricyc2.add_pair('cycle')
trishaded = Picture.new('trishaded')
trioutlined = Picture.new('trioutlined')
trishift = Pair.new((-v11n[0]/3).cm,(-v11n[0]/3).cm)
trishaded.add_drawable(Fill.new(tricyc1).translate(trishift.x,trishift.y).colour(0.7,0.7,0.7))
trishaded.add_drawable(Fill.new(tricyc2).translate(trishift.x,trishift.y).colour(0.7,0.7,0.7))
trioutlined.add_drawable(Draw.new(tricyc1).translate(trishift.x,trishift.y).colour(0.6,0.6,0.6))
trioutlined.add_drawable(Draw.new(tricyc2).translate(trishift.x,trishift.y).colour(0.6,0.6,0.6))


tritesselate = Picture.new('tritesselate')
(-5..5).each do |x|
	(-5..5).each do |y|
		v = M*Vector[x,y]
		tritesselate.add_drawable(Draw.new(trioutlined).translate(v[0].cm,v[1].cm))
	end
end

fig1 = Figure.new
fig1.add_drawable(Clip.new(Square.new, pic).scale(7.93.cm))
fig1.add_drawable(Clip.new(Square.new, tritesselate).scale(7.93.cm))
fig1.add_drawable(Draw.new(trishaded))
fig1.add_drawable(Draw.new(pic))
fig1.add_drawable(Draw.new(tritesselate).translate(8.4.cm, 0))
fig1.add_drawable(Draw.new(pic).translate(8.4.cm, 0))
file.add_figure(fig1)


#####################
# Figure with rectangular tesselating region
#####################
rectcyc = Path.new
rectcyc.add_pair(Pair.new(0.407692307692308.cm,0.561538461538462.cm))
rectcyc.add_pair(Pair.new(-0.592307692307692.cm,0.361538461538461.cm))
rectcyc.add_pair(Pair.new(-0.407692307692308.cm,-0.561538461538462.cm))
rectcyc.add_pair(Pair.new(0.592307692307692.cm,-0.361538461538461.cm))
rectcyc.add_pair('cycle')
rectshaded = Picture.new('rectshaded')
rectoutlined = Picture.new('rectoutlined')
rectshaded.add_drawable(Fill.new(rectcyc).scale(matscale).colour(0.7,0.7,0.7))
rectoutlined.add_drawable(Draw.new(rectcyc).scale(matscale))


recttesselate = Picture.new('recttesselate')
(-5..5).each do |x|
	(-5..5).each do |y|
		v = M*Vector[x,y]
		recttesselate.add_drawable(Draw.new(rectoutlined).translate(v[0].cm,v[1].cm))
	end
end

#fig2 = Figure.new
#fig2.add_drawable(Clip.new(Square.new, pic).scale(7.93.cm))
#fig2.add_drawable(Clip.new(Square.new, recttesselate).scale(7.93.cm))
#fig2.add_drawable(Draw.new(rectshaded))
#fig2.add_drawable(Draw.new(pic))
#fig2.add_drawable(Draw.new(recttesselate).translate(8.4.cm, 0))
#fig2.add_drawable(Draw.new(pic).translate(8.4.cm, 0))
#file.add_figure(fig2)

cliprect = '(-4cm,-3cm)--(4cm,-3cm)--(4cm,3cm)--(-4cm,3cm)--cycle'

fig2 = Figure.new
fig2.add_drawable(Draw.new(rectshaded))
fig2.add_drawable(Clip.new(cliprect, pic))
fig2.add_drawable(Clip.new(cliprect, recttesselate))
fig2.add_drawable(Draw.new(pic))
fig2.add_drawable(Draw.new(recttesselate).translate(0, 0))
fig2.add_drawable(Draw.new(pic).translate(0, 0))
file.add_figure(fig2)



file.compile

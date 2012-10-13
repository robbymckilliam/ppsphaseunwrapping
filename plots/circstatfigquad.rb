require 'rubypost'

require 'matrix'
require 'mathn'

include RubyPost

file = RubyPost::File.new('circstatfigquad')

#def dround(d) return (10000.0*d).round/10000.0 end
def dround(d) return d end

def drawpoly(e, s, func)
	polypath = Path.new
	t = 0
	while(t < e) do
		polypath.add_pair(Pair.new(dround(t).cm, dround(func.call(t)).cm))
		t+=s
	end
	return polypath
end

def drawsamplepoly(e, s, func, picture)
	t = 0
	while(t < e) do
		picture.add_drawable(Draw.new(Circle.new()).scale(0.15.cm).translate(dround(t).cm, dround(func.call(t)).cm))
		t+=s
	end
	return picture
end

def drawcircpoly(e, s, func)
	polypath = Path.new.curved
	t = 0
	while(t < e) do
		a = 2*Math::PI*func.call(t)
		x = dround(t*Math.cos(a)/2.0)
		y= dround(t*Math.sin(a)/2.0)		
		polypath.add_pair(Pair.new(x.cm, y.cm))
		t+=s
	end
	return polypath
end

def drawsamplecircpoly(e, s, func, picture)
	t = 0
	while(t < e) do
    a = 2*Math::PI*func.call(t)
		x = dround(t*Math.cos(a)/2.0)
		y= dround(t*Math.sin(a)/2.0)		
    picture.add_drawable(Draw.new(Circle.new()).scale(0.12.cm).translate(x.cm, y.cm))
		t+=s
	end
	return picture
end

poly = lambda{ |t| 0.4*t**2 - 1.5*t+ 1.5 } 
poly2 = lambda{ |t| -0.1*t**2 - 0*t+ 2.5 } 

polypic = Picture.new('polypic')
polypic.add_drawable('drawdblarrow (5.4cm,0)--(0,0)--(0,4.8cm);')
polypic.add_drawable('label.rt(btex $y$ etex, (0,4.8cm));')
polypic.add_drawable('label.top(btex $t$ etex, (5.4cm,0));')
polypic.add_drawable('label.lft(btex $0.5$ etex, (0,0.5cm));')
polypic.add_drawable('draw (-0.05cm, 0.5cm)--(0.05cm,0.5cm);')
polypic.add_drawable(Draw.new(drawpoly(5.01,0.01, lambda { |t| poly.call(t) })))
polypic.add_drawable(Draw.new(drawpoly(5.01,0.01, lambda { |t| poly2.call(t) })).add_option('dashed evenly'))
drawsamplepoly(5.01,1, lambda { |t| poly.call(t) }, polypic)
drawsamplepoly(5.01,1, lambda { |t| poly2.call(t) }, polypic)

polycircpic = Picture.new('polycircpic')
polycircpic.add_drawable('drawdblarrow (0,-2.7cm)--(0,2.7cm);')
polycircpic.add_drawable('drawarrow (0,0)--(5.4cm,0);')
polycircpic.add_drawable('label.rt(btex $\langle y \rangle$ etex, (0,2.7cm));')
polycircpic.add_drawable('label.top(btex $t$ etex, (5.4cm,0));')
polycircpic.add_drawable('label.lft(btex $0.5$ etex, (0,4*0.5cm));')
polycircpic.add_drawable('draw (-0.05cm, 4*0.5cm)--(0.05cm,4*0.5cm);')
polycircpic.add_drawable('label.lft(btex $-0.5$ etex, (0,-4*0.5cm));')
polycircpic.add_drawable('draw (-0.05cm, -4*0.5cm)--(0.05cm,-4*0.5cm);')
polycircpic.add_drawable(Draw.new(drawpoly(5.02,0.005, lambda { |t| 4*(poly.call(t) - (poly.call(t)).round) })))
polycircpic.add_drawable(Draw.new(drawpoly(5.02,0.005, lambda { |t| 4*(poly2.call(t) - (poly2.call(t)).round) })).add_option('dashed evenly'))
drawsamplepoly(5.02,1, lambda { |t| 4*(poly.call(t) - (poly.call(t)).round) }, polycircpic)
drawsamplepoly(5.02,1, lambda { |t| 4*(poly2.call(t) - (poly2.call(t)).round) } , polycircpic)

spiralpic = Picture.new('spiralpic')
spiralpic.add_drawable('drawarrow (0,0)--(3cm,0);')
spiralpic.add_drawable('label.top(btex $t$ etex, (3cm,0));')
spiralpic.add_drawable('drawarrow (2.7cm,0)..2.7cm*(cosd(20), sind(20))..2.7cm*(cosd(40), sind(40));')
spiralpic.add_drawable('label.urt(btex $\langle y \rangle$ etex, 2.7cm*(cosd(40), sind(40)));')
spiralpic.add_drawable(Draw.new(drawcircpoly(5.02,0.01, lambda { |t| poly.call(t) })))
spiralpic.add_drawable(Draw.new(drawcircpoly(5.02,0.01, lambda { |t| poly2.call(t) })).add_option('dashed evenly'))
drawsamplecircpoly(5.02,1, lambda { |t| poly.call(t) }, spiralpic)
drawsamplecircpoly(5.02,1, lambda { |t| poly2.call(t) } , spiralpic)

fig1 = Figure.new
fig1.add_drawable(Draw.new(polypic))
fig1.add_drawable(Draw.new(polycircpic).translate((6.5).cm,2.5.cm))
#fig1.add_drawable(Draw.new(spiralpic).translate(2.2.cm,(-9.3).cm))
file.add_figure(fig1)


file.compile

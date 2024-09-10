load drifterbetty.mat
load drifterulysses.mat
load drifterinti.mat
load drifterisis.mat
load driftermagellan.mat
load drifternansen.mat
%%
load coastlines
uif = uifigure;
g = geoglobe(uif);
hold(g,'on')
p1 = geoplot3(g,drifterulysses.lat,drifterulysses.lon,[],"LineWidth",2)
p2 = geoplot3(g,drifterbetty.lat,drifterbetty.lon,[],"LineWidth",2)
p3 = geoplot3(g,drifterulysses.lat,drifterulysses.lon,[],"LineWidth",2)
p4 = geoplot3(g,drifterinti.lat,drifterinti.lon,[],"LineWidth",2)
p5 = geoplot3(g,drifterisis.lat,drifterisis.lon,[],"LineWidth",2)
p6 = geoplot3(g,driftermagellan.lat,driftermagellan.lon,[],"LineWidth",2)
p7 = geoplot3(g,drifternansen.lat,drifternansen.lon,[],"LineWidth",2)
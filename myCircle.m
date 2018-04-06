% Draw a circle at (x0,y0) with radius r
function H = myCircle(r,x0,y0,NOP,color)

phi = [0:2*pi/NOP:2*pi];

for ii = 1:NOP+1
    xcircle(ii) = r*cos(phi(ii)) + x0;
    ycircle(ii) = r*sin(phi(ii)) + y0;
end

H = plot(xcircle, ycircle,color);
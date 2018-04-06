function H=circle(center,radius,NOP,linewidth,Color)
%---------------------------------------------------------------------------------------------
% H=CIRCLE(CENTER,RADIUS,NOP,STYLE,LineWidth,Color)
% NOP: Number of points
%----------------------------------------------------

THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
plot(X,Y,Color,'LineWidth',linewidth);
hold on
scatter(center(1),center(2),Color,'filled');
axis square;

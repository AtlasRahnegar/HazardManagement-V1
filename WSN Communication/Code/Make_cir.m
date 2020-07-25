 function [X,Y] = Make_cir(Xcent,Ycent,R)

X = R*cos(linspace(0, 2*pi, 30)) + Xcent;

Y = R*sin(linspace(0, 2*pi, 30)) + Ycent;
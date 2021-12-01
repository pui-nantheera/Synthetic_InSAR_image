function [f]=expcos(x,a,b,c)
%expcos(x,a,b,c)
% 
% calculates a*exp(-bx)*cos(cx)

f = a*exp(-b*x).*cos(c*x);



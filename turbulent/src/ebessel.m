function [f]=ebessel(x,a,r,w)
%ebessel(x,a,r,w)
% 
% calculates a*exp(-x/r)*J0(2*pi*x/w)

f = a*exp(-x/r).*besselj(0,2*pi*x/w);

function [f] = pendiffebessel(bessmod,cvdav)
% function pendiffebessel(bessmod,cvdav)
% 
% returns L2 norm for difference between column 2 of cvdav
% and function a*exp(-x/r)*J0(2*pi*x/w)  [J0 is bessel function]
%
% mod = [b c]'... a = cvdav for x=0
%
% cvdav(:,1) = distance x
% cvdav(:,2) = covariance 
a =cvdav(1,2);
r =bessmod(1);
w =bessmod(2);

f = norm(cvdav(:,2) - ebessel(cvdav(:,1),a,r,w));
%clf
%plot(cvdav(:,1),cvdav(:,2),'g')
%hold on
%plot(cvdav(:,1),expcos(cvdav(:,1),a,b,c),'b')

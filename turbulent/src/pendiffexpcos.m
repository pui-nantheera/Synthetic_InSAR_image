function [f] = pendiffexpcos(abcmod,cvdav)
% function pendiffexpcos(abcmod,cvdav)
% 
% returns L2 norm for difference between column 2 of cvdav
% and function a*exp(-bx)*cos(cx)
%
% mod = [b c]'... a = cvdav for x=0
%
% cvdav(:,1) = distance x
% cvdav(:,2) = covariance 
a =cvdav(1,2);
b =abcmod(1);
c =abcmod(2);

f = norm(cvdav(:,2) - expcos(cvdav(:,1),a,b,c));
%clf
%plot(cvdav(:,1),cvdav(:,2),'g')
%hold on
%plot(cvdav(:,1),expcos(cvdav(:,1),a,b,c),'b')

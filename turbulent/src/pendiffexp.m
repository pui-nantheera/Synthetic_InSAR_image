function [f] = pendiffexp(alphamod,cvdav)
% function pendiffexp(alphamod,cvdav)
% 
% returns L2 norm for difference between column 2 of cvdav
% and function maxvar*exp(-alpha*x)
%
% mod = [alpha]... maxvar = cvdav for x=0
%
% cvdav(:,1) = distance x
% cvdav(:,2) = covariance 
a =cvdav(1,2); %maxvar

f = norm(cvdav(:,2) - a*exp(-alphamod*cvdav(:,1)));
%clf
%plot(cvdav(:,1),cvdav(:,2),'g')
%hold on
%plot(cvdav(:,1),expcos(cvdav(:,1),a,b,c),'b')

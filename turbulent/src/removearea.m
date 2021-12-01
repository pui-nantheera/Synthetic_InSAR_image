function [msk,BW] = removearea(a_old);
%   removearea     - masks out region selected with mouse (assign 0s) 
%
%       run:    a=removearea(a);
%              
%       in:    a_old    old unwrap mask (ones, NaN for not-unwrapped)
%       out:   msk      mask
%
% Function to remove manually an area with unwrap errors.  
% polygon on the image using the mouse.
% input the old unwrap masks  and produces a new one

if nargin==0,help removearea;return;end;


fprintf(1,'Select the polygon using the mouse...\n');

   [BW,Xi,Yi] = roipoly;  %Xi and Yi are the polyn. vortices
% convert binary image to indexed imag
   [BW,map]=gray2ind(BW);
% convert indexed image to intensity image
   BW = ind2gray(BW,map);
% create the mask
   BW(find(BW==0))=1;
   BW(find(BW==255))=0;
% mask out the zero-values
   BW=double(BW);
   msk=BW;
   

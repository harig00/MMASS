function [dloc1,dloc2]=decay(a,d,an,lx,BE,intv)
% [dloc1,dloc2]=decay(a,d,an,lx,BE,intv)
%
% OUTPUT
%
% dloc1   The mean absolute value in the interval
% dloc2   The median absolute value in the interval
%
% Figures out the decay of the wavelet coefficients
% in order to look for a diagnostic
%
% Last modified by fjsimons-at-alum.mit.edu, May 8th 2003


[tsamp,tskol,tsel]=wtimax(a,d,an,lx,BE);

if intv(1)==-12345 & intv(2)==-12345
  intv=BE;
end

for index=1:length(d)
  % First location estimate: mean of the absolute values
  dloc1(index)=mean(abs(d{index}...
			(tskol{index}>intv(1) & tskol{index}<intv(2))));      
  dloc2(index)=median(abs(d{index}...
			(tskol{index}>intv(1) & tskol{index}<intv(2))));
end

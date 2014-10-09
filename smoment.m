function Mw=smoment(cmtline)
% Mw=smoment(cmtline)
% 
% Calculates seismic moment magnitude from a line of the CMT file or else
% directly from the moment tensor. 
%
% INPUT:
%
% cmtline   A 13-entry vector culled from a cmt file, OR:
%           A 6-entry vector which is already the moment tensor
%
% OUTPUT:
%
% Mw        Scalar seismic moment
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012
%
% EXAMPLE:
%
% From the line
% DUR10.2 EX 26  1.98 0.04  0.44 0.03 -2.42 0.03  0.33 0.02  0.13 0.03 -0.82 0.04
% give it
% [26  1.98 0.04  0.44 0.03 -2.42 0.03  0.33 0.02  0.13 0.03 -0.82 0.04]

% Collect the data
if length(cmtline)==13
  ex=cmtline(1);

  Mrr=10^ex*cmtline(2);
  Mtt=10^ex*cmtline(4);
  Mpp=10^ex*cmtline(6);
  Mrt=10^ex*cmtline(8);
  Mrp=10^ex*cmtline(10);
  Mtp=10^ex*cmtline(12);
else
  Mrr=cmtline(1);
  Mtt=cmtline(2);
  Mpp=cmtline(3);
  Mrt=cmtline(4);
  Mrp=cmtline(5);
  Mtp=cmtline(6);
end

% Perform the calculation
Mw=2/3/log(10)*...
   log(1/sqrt(2)*sqrt(Mrr^2+Mtt^2+Mpp^2+...
		      2*(Mrt^2)+2*(Mrp^2)+2*(Mtp^2)))-10.7;  


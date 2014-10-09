function x=d4boxi(f,levels)
% x=D4BOXI(f,[n1 n2 n3])
%
% OBSOLETE, BUT STILL INSTRUCTIVE. 

% Performs a D4 wavelet transform of the three-dimensional array in x,
% whereby the D4 is used in all three directions (first z, then y, then
% x). This transform is NOT periodized, rather, special filters are used
% for the edges. 
%
% INPUT:
%
% f             The three-dimensional array with the wavelet transform
% [n1 n2 n3]    The number of levels in each direction 
%
% OUTPUT:
%
% x             The size(f) matrix with the inversely transformed coefficients
%
% SEE ALSO: D4BOX, D4BOXSTEP, D4BOXSTEPI, and ANGULARD4WT
%
% Written by Ignace Loris (igloris@vub.ac.be) on 22.06.2009
% Last modified by fjsimons-at-alum.mit.edu, 07/15/2009

% Get the levels
defval('levels',[3 3 3])
n1=levels(1);
n2=levels(2);
n3=levels(3);

% Get the coefficients
cofs=d4boxcof;

% Perform the actual inverse wavelet-transform calculation
x=d4box1i(d4box2i(d4box3i(f,n3,cofs),n2,cofs),n1,cofs);

% Subfunctions in each of the dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=d4box1i(x,level,cof)
% starts from a wavelet decomposition up to level
% level in 1 ... log_2(size(x,1))-1
f=x;
for i=level:-1:1 ; f=d4boxstepi(f,i,1,cof); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=d4box2i(x,level,cof)
% starts from a wavelet decomposition up to level
% level in 1 ... log_2(size(x,2))-1
f=x;
for i=level:-1:1 ; f=d4boxstepi(f,i,2,cof); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=d4box3i(x,level,cof)
% starts from a wavelet decomposition up to level
% level in 1 ... log_2(size(x,3))-1
f=x;
for i=level:-1:1 ; f=d4boxstepi(f,i,3,cof); end


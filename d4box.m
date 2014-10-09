function f=d4box(x,levels)
% f=D4BOX(x,[n1 n2 n3])
%
% OBSOLETE, BUT STILL INSTRUCTIVE. 

% Performs a D4 FORWARD wavelet transform of the three-dimensional array
% in x, whereby the D4 is used in all three directions (first 1, then 2,
% then 3). This transform is NOT periodized, rather, special filters are
% used for the edges. The matrix dimensions must be minimally 16 for n1=1,
% 32 for n1=2, etc.
%
% INPUT:
%
% x             The three-dimensional array 
% [n1 n2 n3]    The number of levels in each direction 
%
% OUTPUT:
%
% f             The size(x) matrix with the transform coefficients
%
% EXAMPLE:
%
% vec=randn(128,32,64); nlev=[2 2 2];
% tic; wav=d4box(vec,nlev); toc
% difer(norm3d(vec)-norm3d(wav),6)
% difer(max(max(max(abs(d4boxi(wav,nlev)-vec)))),6)
% rmpath('/u/fjsimons/PROGRAMS/MFILES/wavelets/')
% addpath('/u/fjsimons/MyPapers/WAVELETS/IgnaceLoris')
% rmpath('/u/fjsimons/MyPapers/WAVELETS/IgnaceLoris')
% addpath('/u/fjsimons/PROGRAMS/MFILES/wavelets/')
%
% SEE ALSO: D4BOXI, D4BOXSTEP, D4BOXSTEPI, and ANGULARD4WT
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

% Perform the actual forward wavelet-transform calculation
f=d4box3(d4box2(d4box1(x,n1,cofs),n2,cofs),n3,cofs);

% Subfunctions in each of the dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=d4box1(x,level,cof)
% up to level
% level in 1 ... log_2(size(x,1))-1
f=x;
for i=1:level; f=d4boxstep(f,i,1,cof); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=d4box2(x,level,cof)
% up to level
% level in 1 ... log_2(size(x,2))-1
f=x;
for i=1:level; f=d4boxstep(f,i,2,cof); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=d4box3(x,level,cof)
% up to level
% level in 1 ... log_2(size(x,3))-1
f=x;
for i=1:level; f=d4boxstep(f,i,3,cof); end


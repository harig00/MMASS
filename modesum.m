function fl=modesum(nfl,nn,el)
% fl=MODESUM(nfl,nn,el)
%
% Performs normal-mode summation efficiently, i.e. given 
% f_{nl}, computes sum_{n} f_{nl}
%
% INPUT:
%
% nfl     An array with a generic function depending on n and l
% nn      A same-size array with the mode/branch/overtone numbers
% el      A same-size array with the spherical harmonic degrees
%
% OUTPUT:
%
% fl       A column vector with row index the degree, starting from 0
%
% Last modified by fjsimons-at-alum.mit.edu, 07/10/2012

% All should be turned into vectors, and be same-sized after that
nfl=nfl(:); nn=nn(:); el=el(:);
difer(length(nfl)-length(nn),[],[],NaN)
difer(length(nfl)-length(el),[],[],NaN)

% Perform the summation and output result
fl=full(sum(sparse(nn+1,el+1,nfl),1))';


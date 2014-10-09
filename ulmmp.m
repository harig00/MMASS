function U=ulmmp(l)
% U=ULMMP(l)
%
% The unitary transformation matrix that maps real to complex harmonics,
% for a single spherical harmonic input degree.
%
% INPUT:
%
% l      Spherical harmonic degree
%
% OUTPUT:
%
% U      The unitary conversion matrix
%
% Note that to convert series of expansion coefficients rather than the
% functions themselves, you need the transpose! Note that the orders are
% going down from +l at (1,1) to -l at (m,n).
%
% See also UMMP, CPX2RSH, RSH2CPX
%
% Last modified by fjsimons-at-alum.mit.edu, 02/21/2010

% Make the two diagonals
nwse=[repmat(-i,1,l) sqrt(2) repmat(1,1,l) ];
swne=[(-1).^(l:-1:1) 0       i*(-1).^(1:1:l)];

U=[diag(nwse)+flipud(diag(swne))]/sqrt(2);

% Check unitarity without message
difer(U'*U-eye(size(U)),[],[],NaN)

function U=ummp(L)
% U=UMMP(L)
%
% The unitary block-transformation matrix that maps real to complex
% harmonics, for all degrees 0:L. 
%
% INPUT:
%
% L      Spherical harmonic bandwidth
%
% OUTPUT:
%
% U      Collection of transformation matrices
%
% Note that to convert series of expansion coefficients rather than the
% functions themselves, you need the transpose!  Note that the orders are
% going down from +l at (1,1) to -l at (m,n).
%
% See also ULMMP, CPX2RSH, RSH2CPX
%
% Last modified by fjsimons-at-alum.mit.edu, 02/21/2010

defval('L',3)

% Initialize
U=repmat(0,(L+1)^2,(L+1)^2);

e=0;
for degree=0:L
  b=e+1;
  e=b+(2*degree+1)-1;
  U(b:e,b:e)=ulmmp(degree);
end


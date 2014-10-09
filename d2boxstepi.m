function x=d2boxstepi(f,level,dim,cofs)
% x=d2boxstepi(f,level,dim,cofs)
%
% Performs one INVERSE iteration (from a certain level) of D2 transform
% along a certain dimension; assuming scaling coefficients are in front,
% followed by wavelet coefficients. This transform is an isometry. It is
% equal to the Haar transform.
%
% INPUT:
%
% f             The three-dimensional array with wavelet and scaling functions
% level         The level we are at now [scalar]
% dim           The index identifying the dimension [scalar]
% cofs          The wavelet and scaling filter coefficients
% 
% OUTPUT:
%
% x             The inverse wavelet transform of f
%
% SEE ALSO: D2BOXSTEP, D2BOXCOF
%
% Inspired by Ignace Loris (igloris@vub.ac.be) on 22.06.2009
% Last modified by fjsimons-at-alum.mit.edu, 08/24/2010

% Figure out dimensions
nall=size(f);
if length(nall)==2
  nall(3)=1;
end

% Initialize output, which you need to always take from the previous step 
x=f;

if level==0
  % Do nothing
  return
end

% To move out of 'level' we expand two sets of k coefficients each into
% 2^(n-level+1) coefficients, where 2^n is the dimension of the data set
k=size(f,dim)/2^level;

% The LF-tap filter length
LF=length(cofs.H0);

% Transform low pass coefficients from previous step
% Regular ones (i.e. excludes first LE and last LE)
% Exclude the same number on the left and right
% Compare with the forward transform but most of all work with what
% follows in leaving the right amount of spots untouched.

for i=1:k
  % Isolate the sets of planes in the right dimension
  % Interior, Lowpass  (Scaling Coefficients)
  fplaneslo=[f(dindeks(  i,dim,nall))'];
  % Interior, Highpass  (Wavelet Coefficients)
  fplaneshi=[f(dindeks(k+i,dim,nall))'];
  % And put the convolutions in the SHIFTED UPSAMPLED spots 
  % See under SYNTHESIS, SN p. 129
  % Watch the delay, which is different in D2BOXSTEPI
  % These indices determine the ranges above!
  x(dindeks(2*i-1,dim,nall))=...
      cofs.F0(1)*fplaneslo+cofs.F1(1)*fplaneshi;
  x(dindeks(2*i+0,dim,nall))=...
      cofs.F0(2)*fplaneslo+cofs.F1(2)*fplaneshi;
end

% All other coefficients (wavelet coeff from previous step) remain

function x=d6boxstepi(f,level,dim,cofs)
% x=d6boxstepi(f,level,dim,cofs)
%
% Performs one INVERSE iteration (from a certain level) of D6 transform
% along a certain dimension; assuming scaling coefficients are in front,
% followed by wavelet coefficients. No periodicity but uses wavelets on
% the interval, i.e. filters are different at edges. This transform is an
% isometry. 
% See Cohen, Daubechies and Vial, ACHA 1993.
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
% SEE ALSO: D6BOXSTEP, D6BOX, D6BOXI, D6BOXCOF
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
% The edge filter length
LE=size(cofs.LLO,2);

% Treat the edges with special filters
% Isolate the indices of the left sets of planes
fleftlo=[f(dindeks(  1,dim,nall))'; ...
	 f(dindeks(  2,dim,nall))'; ...
	 f(dindeks(  3,dim,nall))'; ...
	 f(dindeks(  4,dim,nall))'; ...
	 f(dindeks(  5,dim,nall))'];
flefthi=[f(dindeks(k+1,dim,nall))'; ...
	 f(dindeks(k+2,dim,nall))'; ...
	 f(dindeks(k+3,dim,nall))'; ...
	 f(dindeks(k+4,dim,nall))'; ...
	 f(dindeks(k+5,dim,nall))'];

% And put the convolutions in the right spot
stuff=cofs.LFT'*[fleftlo ; flefthi];
for i=1:LE
  x(dindeks(i,dim,nall))=stuff(i,:);
end

% Transform low pass coefficients from previous step
% Regular ones (i.e. excludes first LE and last LE)
% Exclude the same number on the left and right
% Compare with the forward transform but most of all work with what
% follows in leaving the right amount of spots untouched.
if k-4<5; error('Input signal is not long enough'); end
for i=5:k-4
  % Isolate the sets of planes in the right dimension
  % Interior, Lowpass  (Scaling Coefficients)
  fplaneslo=[f(dindeks(  i+1,dim,nall))'; ...
	     f(dindeks(  i+0,dim,nall))'; ...
	     f(dindeks(  i-1,dim,nall))'];
  % Interior, Highpass  (Wavelet Coefficients)
  fplaneshi=[f(dindeks(k+i+1,dim,nall))'; ...
	     f(dindeks(k+i+0,dim,nall))'; ...
	     f(dindeks(k+i-1,dim,nall))'];
  % And put the convolutions in the SHIFTED UPSAMPLED spots 
  % See under SYNTHESIS, SN p. 129
  % Watch the delay, which is different in D4BOXSTEPI
  % These indices determine the ranges above!
  x(dindeks(2*i-1,dim,nall))=...
      cofs.F0([1 3 5])*fplaneslo+cofs.F1([1 3 5])*fplaneshi;
  x(dindeks(2*i+0,dim,nall))=...
      cofs.F0([2 4 6])*fplaneslo+cofs.F1([2 4 6])*fplaneshi;
end

% Isolate the indices of the right sets of planes
frightlo=[f(dindeks(k-4,dim,nall))'; ...
	  f(dindeks(k-3,dim,nall))'; ...
	  f(dindeks(k-2,dim,nall))'; ...
	  f(dindeks(k-1,dim,nall))'; ...
	  f(dindeks(k-0,dim,nall))'];
frighthi=[f(dindeks(2*k-4,dim,nall))'; ...
	  f(dindeks(2*k-3,dim,nall))'; ...
	  f(dindeks(2*k-2,dim,nall))'; ...
	  f(dindeks(2*k-1,dim,nall))'; ...
	  f(dindeks(2*k-0,dim,nall))'];

% And put the convolutions in the right spot
stuff=cofs.RGT'*[frightlo ; frighthi];
for i=1:LE
  x(dindeks(2*k-(LE-i),dim,nall))=stuff(i,:);
end

% All other coefficients (wavelet coeff from previous step) remain


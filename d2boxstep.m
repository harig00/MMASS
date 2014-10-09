function varargout=d2boxstep(x,level,dim,cofs)
% f=D2BOXSTEP(x,level,dim,cofs)
%
% Performs one FORWARD iteration (ending at a certain level of
% decomposition) of D2 transform along a certain dimension; putting
% scaling coefficients in front, followed by wavelet coefficients. This
% transform is an isometry. It is equal to the Haar transform.
%
% INPUT:
%
% x             The three-dimensional array, sized as a power of two
% level         The level we end up at [scalar]
% dim           The index identifying the dimension [scalar]
% cofs          The wavelet and scaling filter coefficients [D2BOXCOF]
%
% OUTPUT:
%
% f             The wavelet transform of x, same dimensions as x
%
% EXAMPLE: My forward is the inverse of my inverse
% d2boxstep('demo1')
%
% SEE ALSO: D2BOXSTEPI, D2BOXCOF
%
% Inspired by Ignace Loris (igloris@vub.ac.be) on 22.06.2009
% Last modified by fjsimons-at-alum.mit.edu, 08/24/2010

if ~isstr(x)

  % Initialize output, which you need to always take from the previous step 
  f=x;

  if level==0
    % Do nothing
    varargout={f};
    return
  end

  % Figure out dimensions
  nall=size(x);
  if length(nall)==2
    nall(3)=1;
  end

  % To move into 'level' we split 2^(n-level+1) coefficients into two sets
  % of k coefficients each, where 2^n is the dimension of the data set
  k=size(x,dim)/2^level;

  % The LF-tap filter length
  LF=length(cofs.H0);

  if k<=2^(LF/2); warning('Input signal is not long enough for reconstruction'); end

  % Exclude the same number on the left and right, namely none
  for i=1:k
    % Isolate the DOWNSAMPLED sets of planes in the right dimension
    xinside=[x(dindeks(2*i+0,dim,nall))'; ...
	     x(dindeks(2*i-1,dim,nall))'];
    % And put the convolutions in the right spot
    % See under ANALYSIS, SN p. 123
    % Interior, Lowpass  (Scaling Coefficients)
    f(dindeks(  i,dim,nall))=cofs.H0*xinside;
    % Interior, Highpass  (Wavelet Coefficients)
    f(dindeks(k+i,dim,nall))=cofs.H1*xinside;
  end

  % All other coefficients (wavelet coeff from previous step) remain
  varargout={f};
elseif strcmp(x,'demo1')

  cofs=d2boxcof; dim=ceil(rand*3);
  % Random-sized array must be at least 2^3 long
  n=ceil(rand*7+3);
  % Edge treatment only accurate to level 0 if it is 2^3 long
  level=ceil(rand*(n-3)); 

  % The random data
  if dim==1
    x=rand([2^n 1   1  ]);
  elseif dim==2
    x=rand([1   2^n 1  ]);
  elseif dim==3
    x=rand([1   1   2^n]);
  end

  % Initial output
  disp(sprintf('\n====== D2BOXSTEP versus D2BOXSTEPI ===== \n'))
  disp(sprintf('n = %i ; lev = %i ; dim = %i',n,level,dim))

  % The forward transform
  xf=d2boxstep(x,level,dim,cofs);
  % The inverse of the forward transform
  xfi=d2boxstepi(xf,level,dim,cofs);
  % The reconstruction error, should be zero
  mae1=mean(mean(mean(abs(x-xfi))));

  % Further output
  % Further output
  disp(sprintf('no preconditioning mean(abs(error)) = %8.3e',mae1))
  disp(sprintf('\n======================================= \n'))
end

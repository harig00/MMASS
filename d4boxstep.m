function varargout=d4boxstep(x,level,dim,cofs)
% f=D4BOXSTEP(x,level,dim,cofs)
%
% Performs one FORWARD iteration (ending at a certain level of
% decomposition) of D4 transform along a certain dimension; putting
% scaling coefficients in front, followed by wavelet coefficients. No
% periodicity but uses wavelets on the interval, i.e. filters are
% different at edges. This transform is an isometry.
% See Cohen, Daubechies and Vial, ACHA 1993.
%
% INPUT:
%
% x             The three-dimensional array, sized as a power of two
% level         The level we end up at [scalar]
% dim           The index identifying the dimension [scalar]
% cofs          The wavelet and scaling filter coefficients [D4BOXCOF]
%
% OUTPUT:
%
% f             The wavelet transform of x, same dimensions as x
%
% EXAMPLE: My forward is the inverse of my inverse
% d4boxstep('demo1')
%
% EXAMPLE: Ignace's forward is the inverse of his inverse
% d4boxstep('demo2')
%
% EXAMPLE: test against Ignace's independent version, 1-D, 1-level, yes/no precon
% d4boxstep('demo3')
%
% SEE ALSO: D4BOX, D4BOXI, D4BOXSTEPI, D4BOXCOF
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

  % Treat the edges with special filters
  % Isolate the indices of the left sets of planes
  xleft=[x(dindeks(1,dim,nall))'; ...
	 x(dindeks(2,dim,nall))'; ...
	 x(dindeks(3,dim,nall))'; ...
	 x(dindeks(4,dim,nall))'; ...
	 x(dindeks(5,dim,nall))'];

  % Left, Lowpass (Scaling Coefficients)
  stuff=cofs.LLO*xleft;
  for i=1:LF/2
    f(dindeks(  i,dim,nall))=stuff(i,:);
  end

  % Left, Highpass (Wavelet Coefficients)
  stuff=cofs.LHI*xleft;
  for i=1:LF/2
    f(dindeks(k+i,dim,nall))=stuff(i,:);
  end

  % See the combination cofs.LFT for D4BOXSTEPI

  % Exclude the same number on the left and right, namely LF/2
  for i=LF/2+1:k-(LF/2-1)
    % Isolate the DOWNSAMPLED sets of planes in the right dimension
    xinside=[x(dindeks(2*i+1,dim,nall))'; ...
	     x(dindeks(2*i+0,dim,nall))'; ...
	     x(dindeks(2*i-1,dim,nall))'; ...
	     x(dindeks(2*i-2,dim,nall))'];
    % And put the convolutions in the right spot
    % See under ANALYSIS, SN p. 123
    % Interior, Lowpass  (Scaling Coefficients)
    f(dindeks(  i,dim,nall))=cofs.H0*xinside;
    % Interior, Highpass  (Wavelet Coefficients)
    f(dindeks(k+i,dim,nall))=cofs.H1*xinside;
  end

  % See the combination cofs.RGT for D4BOXSTEPI

  % Isolate the indices of the right sets of planes
  xright=[x(dindeks(2*k-4,dim,nall))'; ...
	  x(dindeks(2*k-3,dim,nall))'; ...
	  x(dindeks(2*k-2,dim,nall))'; ...
	  x(dindeks(2*k-1,dim,nall))'; ...
	  x(dindeks(2*k-0,dim,nall))'];
  
  % Right, Lowpass (Scaling Coefficients)
  stuff=cofs.RLO*xright;
  for i=1:LF/2
    f(dindeks(  k-(LF/2-i),dim,nall))=stuff(i,:);
  end

  % Right, Highpass (Wavelet Coefficients)
  stuff=cofs.RHI*xright;
  for i=1:LF/2
    f(dindeks(2*k-(LF/2-i),dim,nall))=stuff(i,:);
  end

  % All other coefficients (wavelet coeff from previous step) remain
  varargout={f};
elseif strcmp(x,'demo1')

  cofs=d4boxcof; dim=ceil(rand*3);
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
  disp(sprintf('\n====== D4BOXSTEP versus D4BOXSTEPI ===== \n'))
  disp(sprintf('n = %i ; lev = %i ; dim = %i',n,level,dim))

  % The forward transform
  xf=d4boxstep(x,level,dim,cofs);
  % The inverse of the forward transform
  xfi=d4boxstepi(xf,level,dim,cofs);
  % The reconstruction error, should be zero
  mae1=mean(mean(mean(abs(x-xfi))));

  % The forward preconditioned transform
  xfp=d4boxstep(preconD4(x,[1 1 1],'forward',cofs),level,dim,cofs);
  % The inverse preconditioned of the forward preconditioned transform
  xfpip=preconD4(d4boxstepi(xfp,level,dim,cofs),[1 1 1],'inverse',cofs);
  % The reconstruction error, should be zero
  mae2=mean(mean(mean(abs(x-xfpip))));
  
  % The transpose preconditioned transform
  xfp=d4boxstep(preconD4(x,[1 1 1],'transpose',cofs),level,dim,cofs);
  % The inverse preconditioned of the forward preconditioned transform
  xfpip=preconD4(d4boxstepi(xfp,level,dim,cofs),[1 1 1],'inversetranspose',cofs);
  % The reconstruction error, should be zero
  mae3=mean(mean(mean(abs(x-xfpip))));

  % Further output
  disp(sprintf('no preconditioning mean(abs(error)) = %8.3e',mae1))
  disp(sprintf('precond fwd/inv    mean(abs(error)) = %8.3e',mae2))
  disp(sprintf('precond transpose  mean(abs(error)) = %8.3e',mae3))
  disp(sprintf('\n======================================= \n'))
elseif strcmp(x,'demo2')

  cofs=d4boxcof; dim=ceil(rand*3);
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
  disp(sprintf('\n====== D4INTERVAL versus D4INTERVAL ===== \n'))
  disp(sprintf('n = %i ; lev = %i ; dim = %i',n,level,dim))

  % The forward transform
  xf=d4interval(x,level,false,'forward');
  % The inverse of the forward transform
  xfi=d4interval(xf,level,false,'inverse');
  % The reconstruction error, should be zero
  mae1=mean(mean(mean(abs(x-xfi))));

  % The forward preconditioned transform
  xfp=d4interval(x,level,true,'forward');
  % The inverse preconditioned of the forward preconditioned transform
  xfpip=d4interval(xfp,level,true,'inverse');
  % The reconstruction error, should be zero
  mae2=mean(mean(mean(abs(x-xfpip))));

  % The transpose preconditioned transform
  xfp=d4interval(x,level,true,'transpose');
  % The inverse preconditioned of the forward preconditioned transform
  xfpip=d4interval(xfp,level,true,'inversetranspose');
  % The reconstruction error, should be zero
  mae3=mean(mean(mean(abs(x-xfpip))));

  % Further output
  disp(sprintf('no preconditioning mean(abs(error)) = %8.3e',mae1))
  disp(sprintf('precond fwd/inv    mean(abs(error)) = %8.3e',mae2))
  disp(sprintf('precond transpose  mean(abs(error)) = %8.3e',mae3))
  disp(sprintf('\n======================================= \n'))
elseif strcmp(x,'demo3')
  cofs=d4boxcof; n=ceil(rand*6+3); tolex=9;
  x=rand([2^n 1 1]); level=1; dim=1;

  % Initial output
  disp(sprintf('\n====== D4BOXSTEP versus D4INTERVAL ===== \n'))
  disp(sprintf('n = %i ; lev = %i ; dim = %i',n,level,dim))

  fjs=d4boxstep(x,level,dim,cofs);
  fjsp=d4boxstep(preconD4(x,[1 0],'forward',cofs),level,dim,cofs);
  il=d4interval(x,level,false,'forward');
  ilp=d4interval(x,level,true,'forward');
  difer(fjs-il,tolex,[],sprintf('Agreed to E-%i for size %i',tolex,2^n))
  difer(fjsp-ilp,tolex,[],sprintf('Agreed to E-%i for size %i',tolex,2^n))

  fjsi=d4boxstepi(x,level,dim,cofs);
  ili=d4interval(x,level,false,'inverse');
  fjsip=preconD4(d4boxstepi(x,level,dim,cofs),[1 0],'inverse',cofs);
  ilip=d4interval(x,level,true,'inverse');

  difer(fjsi-ili,tolex,[],sprintf('Agreed to E-%i for size %i',tolex,2^n))
  difer(fjsip-ilip,tolex,[],sprintf('Agreed to E-%i for size %i',tolex,2^n))
  disp(sprintf('\n======================================= \n'))
end


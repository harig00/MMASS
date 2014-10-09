function f=angularD4WT(x,levels,precon,tipe,mr)
% f=ANGULARD4WT(x,[n1 n2],precon,tipe,mr)
%
% Perform D4 (perhaps preconditioned) interval wavelet transform in the
% first two indices of a three-dimensional array up to [n1 n2] scales. 
%
% INPUT:
%
% x         The three-dimensional array, dimensions must be powers of two
% n1,n2     The number of levels in the first two directions
%           n1=0, n2=0 is the identity if precon=[0 0]
% precon    Array of length 2 identifying preconditioning   
%           1 Precondition (NOT an orthogonal transform) so a sequence of
%             ones and a linearly increasing sequence are mapped to a
%             sequence of zeroes in the wavelet bands 
%           0 Don't precondition (YES, this is an orthogonal transform)
%             so a sequence of ones and a linearly increasing sequence are
%             not mapped to all zeroes at the edges in the wavelet bands,
%             i.e. the first and last two coefficients are not zero
% tipe      'forward'|'inverse'|'transpose'|'inversetranspose', always of
%            the same type in both indices
% mr        0 The "original" improper, dimensionally sequential transform
%           1 The "proper" multiresolution transform [default]
% 
% OUTPUT:
%
% f         The three-dimensional array with wavelets and scaling
%           functions. Put scaling coefficients in front, followed by
%           wavelet coefficients. 
%
% SEE ALSO: PRECOND4, D4BOXCOF, D4BOXSTEP, D4BOXSTEPI
% 
% EXAMPLE:
%
% angularD4WT('demo1')
% angularD4WT('demo2')
% angularD4WT('demo3')
%
% Inspired by Ignace Loris (igloris@vub.ac.be) on 26/06/2009
% Last modified by fjsimons-at-alum.mit.edu, 10/30/2010

if ~isstr(x)
  % Get the levels of the decomposition
  defval('levels',[3 3])
  defval('tipe','forward')
  defval('mr',1)

  % Get and check the number of decompositions
  n1=levels(1);
  n2=levels(2);
  if mr==1 && n1~=n2
    error('Number of scales must be equal in both dimensions')
  else
    n=n1;
  end

  % Get the coefficients
  cofs=d4boxcof;

  % Precondition, maybe
  switch tipe
   case {'forward','inversetranspose'}
    x=preconD4(x,precon,tipe,cofs);
  end

  % Prepare for output
  f=x;

  % Now do loop over levels
  switch tipe 
   case {'forward','inversetranspose'}
    if mr==0
      % In the first dimension
      for level=1:n1; f=d4boxstep(f,level,1,cofs); end
      % In the second dimension
      for level=1:n2; f=d4boxstep(f,level,2,cofs); end
    elseif mr==1
      for level=1:n
	% Work on shorter and shorter matrices
	lo1=1:size(f,1)/2^(level-1);
	lo2=1:size(f,2)/2^(level-1);
	% Which makes this look like it's always at the first level
	f(lo1,lo2,:)=d4boxstep(f(lo1,lo2,:),1,1,cofs);
	f(lo1,lo2,:)=d4boxstep(f(lo1,lo2,:),1,2,cofs);
      end
    end
   case {'inverse','transpose'}
    if mr==0
      % In the second dimension
      for level=n2:-1:1; f=d4boxstepi(f,level,2,cofs); end
      for level=n1:-1:1; f=d4boxstepi(f,level,1,cofs); end
    elseif mr==1
      for level=n:-1:1
	% Work on shorter and shorter matrices
	lo1=1:size(f,1)/2^(level-1);
	lo2=1:size(f,2)/2^(level-1);
	% Which makes this look like it's always at the first level
	f(lo1,lo2,:)=d4boxstepi(f(lo1,lo2,:),1,2,cofs);
	f(lo1,lo2,:)=d4boxstepi(f(lo1,lo2,:),1,1,cofs);
      end
    end
  end 

  % Precondition, maybe
  switch tipe
   case {'inverse','transpose'}
    f=preconD4(f,precon,tipe,cofs);
  end
elseif strcmp(x,'demo1')
  % Test Ignace's code for the non-multiresolution case
  whereitsat=fullfile(getenv('MFILES'),'ignaceloris');
  whereiamat=fullfile(getenv('MFILES'),'wavelets');
  x=rand(128,256,32); precon=[1 0];
  ff=angularD4WT(x,[3 4],precon,'forward',0);
  fi=angularD4WT(x,[3 4],precon,'inverse',0);
  fit=angularD4WT(x,[3 4],precon,'inversetranspose',0);
  ft=angularD4WT(x,[3 4],precon,'transpose',0);
  rmpath(whereiamat)
  addpath(whereitsat)
  ffil=angularD4WT(x,[3 4],precon,'forward');
  fiil=angularD4WT(x,[3 4],precon,'inverse');
  fitil=angularD4WT(x,[3 4],precon,'inversetranspose');
  ftil=angularD4WT(x,[3 4],precon,'transpose');
  minmax(ff(:)-ffil(:))
  minmax(fi(:)-fiil(:))
  minmax(fit(:)-fitil(:))
  minmax(ft(:)-ftil(:))
  rmpath(whereitsat)
  addpath(whereiamat)
elseif strcmp(x,'demo2')
  % Test Ignace's code for the multiresolution case
  whereitsat=fullfile(getenv('MFILES'),'ignaceloris');
  whereiamat=fullfile(getenv('MFILES'),'wavelets');
  x=rand(256,256,32); precon=[1 0];
  ff=angularD4WT(x,[4 4],precon,'forward',1);
  fi=angularD4WT(x,[4 4],precon,'inverse',1);
  fit=angularD4WT(x,[4 4],precon,'inversetranspose',1);
  ft=angularD4WT(x,[4 4],precon,'transpose',1);
  rmpath(whereiamat)
  addpath(whereitsat)
  ffil=angularD4MRWT(x,[4 4],precon,'forward');
  fiil=angularD4MRWT(x,[4 4],precon,'inverse');
  fitil=angularD4MRWT(x,[4 4],precon,'inversetranspose');
  ftil=angularD4MRWT(x,[4 4],precon,'transpose');
  minmax(ff(:)-ffil(:))
  minmax(fi(:)-fiil(:))
  minmax(fit(:)-fitil(:))
  minmax(ft(:)-ftil(:))
  rmpath(whereitsat)
  addpath(whereiamat)
elseif strcmp(x,'demo3')
  % Looks like I hadn't checked before that we get perfect reconstruction
  x=rand(256,256,32); precon=[1 1];
  disp(sprintf('norm of   x is %8.3f',norm(x(:))))
  ff=angularD4WT(x,[4 4],precon,'forward',1);
  disp(sprintf('norm of  ff is %8.3f',norm(ff(:))))
  xpr=angularD4WT(ff,[4 4],precon,'inverse',1);
  disp(sprintf('norm of xpr is %8.3f',norm(xpr(:))))
  difer(xpr-x,6)
end


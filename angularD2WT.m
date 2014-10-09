function f=angularD2WT(x,levels,tipe,mr)
% f=ANGULARD2WT(x,[n1 n2],tipe,mr)
%
% Perform D2 wavelet transform in the first two indices of a
% three-dimensional array up to [n1 n2] scales. 
%
% INPUT:
%
% x         The three-dimensional array, dimensions must be powers of two
% n1,n2     The number of levels in the first two directions
%           n1=0, n2=0
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
% Inspired by Ignace Loris (igloris@vub.ac.be) on 26/06/2009
% Last modified by fjsimons-at-alum.mit.edu, 11/03/2010

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
  cofs=d2boxcof;

  % Prepare for output
  f=x;

  % Now do loop over levels
  switch tipe 
   case {'forward','inversetranspose'}
    if mr==0
      % In the first dimension
      for level=1:n1; f=d2boxstep(f,level,1,cofs); end
      % In the second dimension
      for level=1:n2; f=d2boxstep(f,level,2,cofs); end
    elseif mr==1
      for level=1:n
	% Work on shorter and shorter matrices
	lo1=1:size(f,1)/2^(level-1);
	lo2=1:size(f,2)/2^(level-1);
	% Which makes this look like it's always at the first level
	f(lo1,lo2,:)=d2boxstep(f(lo1,lo2,:),1,1,cofs);
	f(lo1,lo2,:)=d2boxstep(f(lo1,lo2,:),1,2,cofs);
      end
    end
   case {'inverse','transpose'}
    if mr==0
      % In the second dimension
      for level=n2:-1:1; f=d2boxstepi(f,level,2,cofs); end
      for level=n1:-1:1; f=d2boxstepi(f,level,1,cofs); end
    elseif mr==1
      for level=n:-1:1
	% Work on shorter and shorter matrices
	lo1=1:size(f,1)/2^(level-1);
	lo2=1:size(f,2)/2^(level-1);
	% Which makes this look like it's always at the first level
	f(lo1,lo2,:)=d2boxstepi(f(lo1,lo2,:),1,2,cofs);
	f(lo1,lo2,:)=d2boxstepi(f(lo1,lo2,:),1,1,cofs);
      end
    end
  end 
end

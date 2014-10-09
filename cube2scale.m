function varargout=cube2scale(N,levels,mr,meth)
% [vwlev,vwlevs,vwide]=cube2scale(N,[n1 n2],mr,meth)
%
% After a wavelet transform of a single-depth cubed sphere over all
% faces using ANGULARDXWT (where X=2,4,6), separates the channels into
% wavelet and scaling coefficients at different scales
%
% INPUT:
%
% N        The power of the dyadic subdivision for the cubed sphere
%          ... if [N M] can handle rectangular situations
% n1,n2    The number of levels in the first two dimensions
% mr       0 The "original" improper, dimensionally sequential transform
%          1 The "proper" multiresolution transform [default]
% meth     1 A fast and direct level mr=1 calculation method [default]
%          2 A slow but explicit level mr=1 calculation method
%
% OUTPUT:
% 
% vwlev    A single-face map identifying the decomposition level
% vwlevs   A whole-cube map identifying the decomposition level
% vwide    A single-face map of scaling (0) or wavelet (1) coefficients
%          This is yet to come - we haven't needed it yet
%
% EXAMPLE:
%
% vwlev1=cube2scale(8,[3 3],1,1);
% vwlev2=cube2scale(8,[3 3],1,2);
%
% SEE ALSO:
%
% The function DNUMS which is the one-dimensional equivalent, in a way
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/18/2013

defval('N',4)
defval('mr',1)
defval('levels',[3 3])

% Get and check the number of decompositions
n1=levels(1);
if length(levels)>1
  n2=levels(2);
else
  n2=n1;
end
if mr==1 && n1~=n2
  error('Number of scales must be equal in both dimensions')
else
  J=n1;
end

if length(N)==2
  M=N(1);
  N=N(2);
elseif length(N)==1
  [M,N]=deal(N);
end

% Define method
defval('meth',1)

% Now do loop over levels
if mr==1
  switch meth
   case 1
    vwlev1=min(J,[M-ceil(log2(1:2^M))+1]);
    vwlev2=min(J,[N-ceil(log2(1:2^N))+1]);
    [a,b]=meshgrid(vwlev2,vwlev1);
    vwlev=min(a,b);
   case 2
    % Let's go heuristically with Haar to watch the cancellation, etc
    f=ones(2^M,2^N);

    % Prepare for output
    fm=size(f,1);
    fn=size(f,2);

    for level=1:J
      lo1=1:fm/2^(level-1);
      lo2=1:fn/2^(level-1);
      f(lo1,lo2)=[repmat(level-0,fm/2^level,fn/2^(level-1)) ;...
		  repmat(level-1,fm/2^level,fn/2^(level-1))];
      f(lo1,lo2)=f(lo1,lo2)+...
	  [repmat(level+1,fm/2^(level-1),fn/2^level) ...
	   repmat(level+0,fm/2^(level-1),fn/2^level)];
    end
    % The scales are now identified by 
    vwlev=min(J,ceil(f/2));
  end
  % imagesc(vwlev); grid on
  % set(gca,'xtick',2.^[N-J:N]+0.5); set(gca,'ytick',2.^[N-J:N]+0.5)
elseif mr==0
  % Let's go heuristically with Haar to watch the cancellation, etc
  f=ones(2^M,2^N);
  
  % Prepare for output
  fm=size(f,1);
  fn=size(f,2);
  
  % In the first dimension
  for level=1:n1
    lo1=1:fm/2^(level-1);
    lo2=1:fn;
    f(lo1,lo2)=[repmat(level,fm/2^level,fn) ;...
		repmat(level-1,fm/2^level,fn)];
  end
  % In the second dimension
  for level=1:n2
    lo1=1:fm;
    lo2=1:fn/2^(level-1);
    f(lo1,lo2)=f(lo1,lo2)+...
	[repmat(level+1,fm,fn/2^level) ...
	 repmat(level+0,fm,fn/2^level)];
  end
  % The scales are now identified by 
  vwlev=min(J,ceil(f/2));
end

if nargout>=2
  % Make the whole-cube map
  vwlevs=repmat(vwlev,[1 1 6]);
else
  vwlevs=NaN;
end

% Prepare output
varns={vwlev,vwlevs};
varargout=varns(1:nargout);

% Junk that didn't work
% This here is no longer accurate sometimes
% vwlev=sqrt(vwlev(:)*vwlev(:)');
% vwlev=nextpow2(floor(vwlev));
% Got to find a solution for this somehow
% vwlev=floor((a+b)/2);

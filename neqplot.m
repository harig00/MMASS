function varargout=neqplot(glon,glat,mapt,varargin)
% NEQPLOT(glon,glat,mapt)
% [LON,LAT,phandle]=NEQPLOT(glon,glat,mapt)
% [LON,LAT,phandle]=NEQPLOT(glon,glat,mapt,'Property','Value')
%
% For a square, regular grid and a parameterization
% from clutser, plots the associated grid. 
% The mapping matrix just has running indices for each separate domain.
%
% 'glon' and 'glat' are just vectors, before meshgrid
%
% EXAMPLE:
%
% N=64;
% mat=peaks(N);
% mapt=clutser(mat,8,4,2);
% glon=0:N;
% glat=N:-1:0;
% imagef([1/2 N-1/2],[N-1/2 1/2],mat) ; hold on
% neqplot(glon,glat,mapt)
%
% Last modified by fjsimons-at-alum.mit.edu, 09/11/2007

if prod(size(glon))-1 ~= size(mapt,2) ...
      | prod(size(glat))-1 ~= size(mapt,1) 
  error('Matrix and vectors incompatible')
end

[LON,LAT]=meshgrid(glon,glat);

lini=~diff(mapt,1,1);
linj=~diff(mapt,1,2);         

linj=[zeros(size(linj(:,1))) linj];
lini=[zeros(size(lini(1,:))) ; lini];

[deli,delj]=find(lini.*linj);

rind=(delj-1)*size(LON,1)+deli;

LON(rind)=NaN;
LAT(rind)=NaN;

if ~nargout | nargout==3
  if nargin==3
    phan=fridplot(LON,LAT);
  else
        phan=fridplot(LON,LAT,varargin{1},varargin{2});
  end
end

strout={ 'LON','LAT','phan'};
for index=1:nargout
  eval([ 'varargout{index}=',strout{index},';'])
end

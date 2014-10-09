function varargout=cube2lola(v,ins)
% [lons,lats,vd,lon,lat]=cube2lola(v,ins)
%
% Transforms a cubed-sphere model to geographical coordinates by
% fine-scale interpolation to a regular grid
%
% INPUT:
%
% v            NxN(x6) matrix with model values
% ins          If v is not NxNx6 then it's this index
%              or if v is NxNx6 then it only takes this index
%
% OUTPUT:
%
% lons,lats    The longitude-latitude grid [radians]
% vd           Corresponding data values
% lon,lat      The longitude-latitude values of the givens [radians]
%
% Last modified by fjsimons-at-alum.mit.edu, 06/28/2011

defval('ins',1:6)

% Find the coordinates of the cubed sphere that are being used
[x,y,z]=cube2sphere(nextpow2(size(v,1)));

% Turn all them into regular spherical coordinates
[lon,lat]=cart2sph(x,y,z);
lon(lon<0)=lon(lon<0)+2*pi;

% Now make target regular geographical coordinates
lons=linspace(0   ,2*pi  ,4*size(v,2));
lats=linspace(pi/2, -pi/2,2*size(v,1));
[lons,lats]=meshgrid(lons,lats);

% If only one face
if length(ins)==1
  lon=lon(:,:,ins);
  lat=lat(:,:,ins);
  difer(prod(size(v))-prod(size(lon)),[],[],NaN)
  difer(prod(size(v))-prod(size(lat)),[],[],NaN)
  % Remember that this was flipped, see PLOTONCUBE
  v=flipud(v);
else
  % Remember that this was flipped, see PLOTONCUBE
  for i=ins; v(:,:,i)=flipud(v(:,:,i)); end
end

warning off MATLAB:griddata:DuplicateDataPoints 
% Perform the interpolation
vd=griddata(lon(:),lat(:),v(:),lons,lats);
warning on MATLAB:griddata:DuplicateDataPoints 

% Output
varns={lons,lats,vd,lon,lat};
varargout=varns(1:nargout);


function [lonm,latm,vd]=cube2moll(v)
% [lonm,latm,vd]=cube2moll(v)
%
% Transforms a cubed-sphere model to Mollweide projection by
% fine-scale interpolation to a regular grid
%
% INPUT:
%
% v            NxNx6 matrix with model values
%
% OUTPUT:
%
% lonm,latm    Mollweide longitudes and latitudes
% vd           Corresponding data values
%
% Last modified by fjsimons-at-alum.mit.edu, 06/28/2010

% Find the coordinates of the cubed sphere that are being used
[x,y,z]=cube2sphere(nextpow2(size(v,1)));
% Turn them into regular spherical coordinates
[lon,lat]=cart2sph(x,y,z);
lon(lon<0)=lon(lon<0)+2*pi;
% Make Mollweide coordinates fit for the purpose
lons=linspace(0,2*pi,4*size(v,2));
lats=linspace(pi/2,-pi/2,4*size(v,1));
[lons,lats]=meshgrid(lons,lats);
[lonm,latm]=mollweide(lons,lats,pi);

% Remember that this was flipped, see PLOTONCUBE
for i=1:6; v(:,:,i)=flipud(v(:,:,i)); end

warning off MATLAB:griddata:DuplicateDataPoints 
% Perform the interpolation ON THE SPHERE!
vd=griddata(lon(:),lat(:),v(:),lons,lats);
warning on MATLAB:griddata:DuplicateDataPoints 

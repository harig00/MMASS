function [names,lon,lat,par]=locident(file,lonlat,tol)
% [names,lon,lat,par]=LOCIDENT(file,lonlat,tol)
%
% For a file containing 
%
% name lat lon parameter
% 
% returns those names for whom the target coordinate pair is within some
% tolerance.
%
% INPUT:
%
% file      The file in question
% lonlat    The target coordinate pair
% tol       The distance tolerance, in km
%
% OUTPUT:
%
% names     Directory names of those close
% lon       Longitudes of those close
% lat       Latitudes of those close
% par       The parameter of those close
%
% Last modified by fjsimons-at-alum.mit.edu, 13.04.2006

defval('lonlat',[242.2243   35.8791])
defval('lonlat',[243.75 34.61])
defval('file','/home/fjsimons/EALARMS/programming/namelocmag');
defval('tol',40)
nlm=load(file);

gcdkm=grcdist(lonlat,[nlm(:,2) nlm(:,3)]);
pos=1:size(nlm);
names=sprintf('%15.6f',nlm(pos(gcdkm<tol),1));
names=reshape(names',15,length(names)/15)';

lon=nlm(pos(gcdkm<tol),2)+360;
lat=nlm(pos(gcdkm<tol),3);
par=nlm(pos(gcdkm<tol),4);



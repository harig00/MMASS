function gridcheck(lonlat,c11,cmn)
% GRIDCHECK(lonlat,c11,cmn)
% 
% Checks if the grid definition is adequate to contain all events and  stations
%
% INPUT:
%
% lonlat    A list of longitudes and latitudes [degrees, longitudes 0->360]
%           [lon1 lat1 lon2 lat2] for stations and events
% c11       [x,y]/[lon,lat]-coordinates of the upper left grid corner [degrees]
% cmn       [x,y]/[lon,lat]-coordinates of the bottom right corner [degrees]
%
% Last modified by fjsimons-at-alum.mit.edu, 06/12/2007

[m,n]=size(lonlat);
if ~all([min(lonlat(:,1),lonlat(:,3)) min(lonlat(:,2),lonlat(:,4))]>...
      repmat([c11(1) cmn(2)],m,1) & ...
      [max(lonlat(:,1),lonlat(:,3)) max(lonlat(:,2),lonlat(:,4))]<...
      repmat([cmn(1) c11(2)],m,1))
  error('Make your grid bigger!')
end


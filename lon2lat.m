function [lonx,latx]=lon2lat(lon1lat1,lon2lat2,lon,fullorsens)
% [lonx,latx]=LON2LAT([lon1 lat1],[lon2 lat2],lon,1)
%        For all given longitudes
% [lonx,latx]=LON2LAT([lon1 lat1],[lon2 lat2],lon,2)
%       Only for those within the anchor points (between max and min longitude)
%
% Given a vector 'lon' of  longitude(s), calculates latitude(s) (in degrees)
% lying on the great circle going through ['lon1' 'lat1'] and ['lon2' 'lat2'].
% Now accepts more than one great circle in column form.
%
% EXAMPLE 1:
%
% % Do Example 2 of lat2lon and then
% [lona,lata]=lon2lat([183 -3.47],[118 -35],lonx,2);
% hold on; plot(lona,lata,'ro')
%
% EXAMPLE 2:
% 
% % Or to illustrate the multiple great-circle function
% [lona,lata]=...
% lon2lat([183 -3.47; 120 -12 ; 130 -5],[123 -35 ; 160 -40 ; 180 -42],lonx,2);
% hold on; plot(lona,lata); 
% plot([183  120 130 123 160 180],[-3.47 -12 -5 -35 -40 -42],'rs')
%
% See also GRCDEMO, LAT2LON
%
% Written by fjsimons-at-alum.mit.edu, 07/22/1999
% Last modified by fjsimons-at-alum.mit.edu, 07/24/2014

% Conversion to radians
lon=lon(:)*pi/180;
lon1lat1=lon1lat1*pi/180;
lon2lat2=lon2lat2*pi/180;

[lon1,lat1]=deal(lon1lat1(:,1),lon1lat1(:,2));
[lon2,lat2]=deal(lon2lat2(:,1),lon2lat2(:,2));

% Number of great circles
n=length(lon1);

% Keep within the anchor points of great-cricle 
% that spans the widest longitude; re-adjust later
if fullorsens==2
  lonx=lon(lon<=max([lon1(:) ; lon2(:)])+2*eps ...
      & lon>=min([lon1(:) ; lon2(:)])+2*eps);
else
  lonx=lon;
end

% Number of crossings
m=length(lonx);

% Make convenient matrix expansions
if n==1; m=1;  end
LAT1=repmat(lat1',m,1); LON1=repmat(lon1',m,1); 
LAT2=repmat(lat2',m,1); LON2=repmat(lon2',m,1);
LONX=repmat(lonx,1,n);

% Do calculation
latx=atan((sin(LAT1).*cos(LAT2).*sin(LONX-LON2)...
     -sin(LAT2).*cos(LAT1).*sin(LONX-LON1))...
     ./(cos(LAT1).*cos(LAT2).*sin(LON1-LON2)));

% Re-adjustment: each great circle within its bounds
for index=1:n
  latx(lonx>max(lon1(index),lon2(index)) ...
      | lonx<min(lon1(index),lon2(index)),index)=NaN;
end

% Conversion to degrees
latx=latx*180/pi;
lonx=lonx*180/pi;

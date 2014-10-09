function [lonx,latx]=lat2lon(lon1lat1,lon2lat2,lat,fullorsens)
% [lonx,latx]=LAT2LON([lon1 lat1],[lon2 lat2],lat,1)
%     A complete great circle
% [lonx,latx]=LAT2LON([lon1 lat1],[lon2 lat2],lat,2)
%     Only within the longitudes of the anchor points
%
% Given a vector of latitude(s), calculates longitude(s) (in degrees)
% lying on the great circle going through [lon1 lat1] and [lon2 lat2].
% Any given great circle (excepting one over the poles) crosses each
% meridian once and only once. However, any given great circle has a maximum
% latitude reached at its apex. It crosses lower latitudes twice and higher
% latitudes never. Output are those longitudes that are
% within the given great circle longitudes. Returns latitudes and longitudes
% at which the path crosses the great circle (in degrees).
% In the case of an incomplete great circle, if either of the two points is
% on the apex, that will not be included; it needs to be done manually
% with 'apex2'. For the full great circle, the apex and its antipode will be
% included.
%
% EXAMPLE 1 (if full great circle)
%
% [lonx,latx]=lat2lon([0 50]*rand,[360 50]*rand,linspace(-90,90,500),1);
% plot(lonx,latx,'b+') ; axis equal 
% axis([0 350 -90 90]); grid;  plotcont([0 90],[360 -90]); hold off; 
%
% EXAMPLE 2 (standard)
%
% [lonx,latx]=lat2lon([183 -3.47],[118 -35],linspace(-90,90,500),2);
% plot(lonx,latx,'b+') 
% grid;  plotcont([0 90],[360 -90]); axis equal; axis([100 200 -50 0]);
%
% See also GRCDEMO
%
% Written by fjsimons-at-alum.mit.edu, 07/22/1999
% Last modified by fjsimons-at-alum.mit.edu, 07/24/2014

% Conversion to radians
lat=lat(:)*pi/180;
lon1lat1=lon1lat1*pi/180;
lon2lat2=lon2lat2*pi/180;
[lon1,lat1]=deal(lon1lat1(:,1),lon1lat1(:,2));
[lon2,lat2]=deal(lon2lat2(:,1),lon2lat2(:,2));

% Number of great circles
n=length(lon1);

% Does our latitude array contain the apex? Then only crosses once.
%[latmax,lonmax]=apex2(lon1lat1*180/pi,lon2lat2*180/pi)
%if any(sum(lat==latmax*pi/180))
%  apexdet=1;
%  disp('Apex detected')
%else
%  apexdet=0;
%end

m=length(lat);

% Make convenient matrix expansions
if n==1; m=1;  end % Unchanged if only one path
LAT1=repmat(lat1',m,1); LON1=repmat(lon1',m,1); 
LAT2=repmat(lat2',m,1); LON2=repmat(lon2',m,1);
LAT=repmat(lat,1,n);

DLON=LON1-LON2;

A=sin(LAT1).*cos(LAT2).*sin(DLON).*cos(LAT);
B=sin(LAT1).*cos(LAT2).*cos(DLON).*cos(LAT)-cos(LAT1).*sin(LAT2).*cos(LAT);
C=cos(LAT1).*cos(LAT2).*sin(LAT).*sin(DLON);

cmat=C./sqrt(A.^2+B.^2);
lonex=atan2(B,A);
cacos=acos(cmat);

% No crossing (should have been intercepted by min/max protection above)
nox=~~imag(cacos);
LAT(nox)=NaN;
lonex(nox)=NaN;
cacos(nox)=NaN;

% First crossing
lon_1=mod(LON1+lonex+cacos+pi,2*pi)-pi;
lon_1=lon_1+(lon_1<0)*2*pi;

% Second crossing
lon_2=mod(LON1+lonex-cacos+pi,2*pi)-pi;
lon_2=lon_2+(lon_2<0)*2*pi;

% Stay within longitude range if not full great circle
% Allow for multiple latitude crossings but only return one vector
if fullorsens==2
  stw1=lon_1<=max([lon1(:) ; lon2(:)])+2*eps & lon_1>=min([lon1(:) ; lon2(:)])-2*eps;
  stw2=lon_2<=max([lon1(:) ; lon2(:)])+2*eps & lon_2>=min([lon1(:) ; lon2(:)])-2*eps;
  %latx=[latx(stw1) ; latx(stw2)];
  LAT(~(stw1 | stw2))=NaN;
  lon_1(~(stw1))=NaN;
  lon_2(~( stw2))=NaN;
  LAT=[LAT ; LAT];
  lon_1=[lon_1 ; lon_2];
else
  [stw1,stw2]=deal([]);
end

% Apex protection
%if fullorsens~=2
%  if apexdet==1
%      disp('Apex properly handled')
%    [latx,ind]=sort([latx -latmax latmax]);
%    lon_1=[lon_1 (lonmax-180)+((lonmax-180)<0)*2*pi lonmax];
%    lon_2=[lon_2 NaN NaN];
%    lon_1=lon_1(ind);
%    lon_2=lon_2(ind);
%  end  
%end

% No crossing whatsoever (new)
if any(stw1 | stw2)
  if n==1
    getrid=~isnan(LAT)&~isnan(lon_1);
    LAT=LAT(getrid);
    lon_1=lon_1(getrid);
    % Sort according to gcdist (to be fancy)      
    [mx,nx]=size(lon_1);
    dist=grcdist([lon_1 LAT]*180/pi,repmat([lon1 lat1]*180/pi,mx,1));
    [dist,sind]=sort(dist);
    lon_1=lon_1(sind);
    LAT=LAT(sind);
  else
    % Readjustment: each great circle within its bounds
    for index=1:n
      LAT(lon_1(:,index)>max(lon1(index),lon2(index))...
	  | lon_1(:,index)<min(lon1(index),lon2(index)),...
	  index)=NaN;
    end
  end
  
  % Conversion to degrees
  latx=LAT*180/pi;
  lonx=lon_1*180/pi;
  
else
  [lonx,latx]=deal([]);
end

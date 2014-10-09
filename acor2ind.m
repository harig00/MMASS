function [celnr,rownr,colnr]=acor2ind(lat,dlon,nmr,c11,lon1lat1)
% [celnr,rownr,colnr]=ACOR2IND(lat,dlon,nmr,c11,[lon1 lat1])
%
% Matches a physical coordinate to a sequential equal-area grid cell index
% 
% lat        latitudes describing the grid, column vector [degrees]
% dlon       longitudes describing the grid, column vector [degrees]
% nmr        the number of cells per row in the grid
% c11        [x,y]/[lon,lat]-coordinates of the upper left grid corner [degrees]
% lon1,lat1  latitudes and longitudes of the points whose indices you want
%
% EXAMPLE:
%
% c11=[110 10]; cmn=[180 -50];
% [lat,dlon,c,nmr]=authalic(c11,cmn,2,5);
% clf; colormap spring
% eqplot(lat,dlon,[c11 cmn]); hold on
% lon1=c11(1)+abs(randn(10,1)*30);
% lat1=c11(2)-abs(rand(10,1)*30);
% [celnr,rownr,colnr]=acor2ind(lat,dlon,nmr,c11,[lon1 lat1]);
% fillauth(lat,dlon,nmr,c11,celnr,rand(10,1))
% hold on; plot(lon1,lat1,'ks')
%
% Last modified by fjsimons-at-alum.mit.edu, 06/13/2007

[lon1,lat1]=deal(lon1lat1(:,1),lon1lat1(:,2));
lat=lat(:);

cumonin=cumsum([0 ; nmr(:)]);

rownr=sum(repmat(lat,1,length(lat1))...
    >repmat(lat1(:)',length(lat),1),1)';

colnr=ceil((lon1(:)-c11(1))./dlon(rownr));

celnr=cumonin(rownr)+colnr;

% Checks that the reverse works, too
[rc,cc]=aind2sub(lat,dlon,nmr,celnr);
if ~all([rc cc]==[rownr colnr]); error ; end

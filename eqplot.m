function varargout=eqplot(lat,dlon,c11cmn,varargin)
% ah=EQPLOT(lat,dlon,[c11 cmn])
%
% Plots an equal-area grid as calculated by AUTHALIC
%
% INPUT:
%
% lat       latitudes describing the grid, column vector [degrees]
% dlon      longitudes describing the grid, column vector [degrees]
% c11       [x,y]/[lon,lat]-coordinates of the upper left grid corner [degrees]
% cmn       [x,y]/[lon,lat]-coordinates of the bottom right corner [degrees]
%
% OUTPUT:
%
% ah        A set of graphics handles:
%              ah(1:length(lat)) % handles representing lines of latitude
%              ah(length(lat)+1:end) % handles for the longitudinal divisions
%
% EXAMPLE:
%
% c11=[110 10]; cmn=[180 -50];
% [lat,dlon,c,nmr]=authalic(c11,cmn,1,5);
% a=eqplot(lat,dlon,[c11 cmn]); set(a(length(lat)+1:end),'color','b')
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/12/2007

[c11,cmn]=deal(c11cmn(1:2),c11cmn(3:4));

% Plots lines of constant latitude and returns handles to it
lats=plot([c11(1) cmn(1)],[lat lat],'k');
hold on

if nargin==3
  lonm=repmat(NaN,2*size(dlon,1),size(c11(1):min(dlon):cmn(1),2));
  for index=1:length(lat)-1
    divs=[c11(1):dlon(index):cmn(1)];
    lonm(2*index-1,1:length(divs))=divs;
    lonm(2*index,1:length(divs))=divs;
  end
  latm=repmat(lat(2:end-1)',2,1);
  latm=latm(:);
  latm=repmat([lat(1) ; latm ; lat(end)],1,size(lonm,2));
  latm(isnan(lonm))=NaN;
  % So a single handle returns a curvy longitudinal line
  hands=plot(lonm,latm,'k');
else
  for index=1:length(lat)-1
    ah=plot(repmat([c11(1):dlon(index):cmn(1)],2,1),...
	[lat(index) lat(index+1)],'k');
    hands=[hands ; ah];
  end
end

axis image
hold off

if nargout~=0
  varargout{1}=[lats ; hands];
end


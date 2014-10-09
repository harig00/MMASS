function varargout=fillauth(lat,dlon,nmr,c11,celnr,colvec,varargin)
% FILLAUTH(lat,dlon,nmr,c11,celnr,colvec)
% FILLAUTH(lat,dlon,nmr,c11,celnr,colvec,caksis)
% fbs=FILLAUTH(lat,dlon,nmr,c11,celnr,colvec,...)
%
% Fills the cells of an equal-area grid with prescribed colors
%
% INPUT:
%
% lat       latitudes describing the grid, column vector [degrees]
% dlon      longitudes describing the grid, column vector [degrees]
% c11       [x,y]/[lon,lat]-coordinates of the upper left grid corner [degrees]
% cmn       [x,y]/[lon,lat]-coordinates of the bottom right corner [degrees]
% nmr       the number of cells per row in the grid
% celnr     a vector with cell numbers that you would like to color in
% colvec    a color string, or a color index into the current color map
% caksis    an optional color range
%
% OUTPUT:
%
% fbsb      graphics handles to the filled patches
%
% EXAMPLES:
% 
% c11=[110 0]; cmn=[180 -50];
% [lat,dlon,c,nmr]=authalic(c11,cmn,1,0.5);
% eqplot(lat,dlon,[c11 cmn]); hold on
%% And then try one of the lines below
% fillauth(lat,dlon,nmr,c11,[1:13 178:211 1324:1473],'b')
% fillauth(lat,dlon,nmr,c11,1:sum(nmr),2)
% fillauth(lat,dlon,nmr,c11,1:sum(nmr),rand(sum(nmr),1))
% fillauth(lat,dlon,nmr,c11,1:sum(nmr),rand(sum(nmr),1)')
% fillauth(lat,dlon,nmr,c11,1:sum(nmr),1:sum(nmr))
%
% Last modified by fjsimons-at-alum.mit.edu, 06/07/2007

colmap=colormap;
celnr=celnr(:);

% What row and column are they on?
[rownr,colnr,top,bot,lef,rig]=aind2sub(lat,dlon,nmr,celnr,c11);

% Color indices
if nargin==7
  colvec=ceil(dat2col(colvec,varargin{1},1));
else
 fbs=fillbox([lef rig top bot],colvec);
end

if nargout
  varargout{1}=fbs;
end


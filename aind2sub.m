function [rownr,colnr,varargout]=aind2sub(lat,dlon,nmr,celnr,varargin)
% [rownr,colnr]=AIND2SUB(lat,dlon,nmr,celnr)
% [rownr,colnr,t,b,l,r]=AIND2SUB(lat,dlon,nmr,celnr,c11)
%
% Returns row and 'column' index given a sequential equal-area grid cell number 
%
% INPUT:
%
% lat         latitudes describing the grid, column vector [degrees]
% dlon        longitudes describing the grid, column vector [degrees]
% nmr         the number of cells per row in the grid
% celnr       the requested running cell index vector
% c11         [x,y]/[lon,lat]-coordinates of the upper left grid corner [degrees]
% 
% OUTPUT:
%
% rownr       row numbers of the requested sequential cell numbers
% colnr       'column' numbers of the requested sequential cell numbers
% [t,b,l,r]   top,bot,right and left BOUNDARIES of the cell if the coordinates;
%             this output only if C11, the BOUNDARIES OF THE first cell are input
%
% EXAMPLE:
%
% c11=[110 0]; cmn=[180 -50];
% [lat,dlon,c,nmr]=authalic(c11,cmn,3,3);
% eqplot(lat,dlon,[c11 cmn]); hold on
%% Fill the 104th cell with a red color
% celnr=104;
% fillauth(lat,dlon,nmr,c11,celnr,'r')
%% And verify the row and the 'column' numbers of this particular cell
% [rn,cn,t,b,l,r]=aind2sub(lat,dlon,nmr,celnr,c11);
% title(sprintf('Cell number %i is row %i and column %i',celnr,rn,cn))
%% Verify the spatial coordinates of this particular cell
% set(gca,'xtick',[l r],'ytick',[b t],'xgrid','on','ygrid','on')
%
% See also ACOR2IND, ASUB2IND
%
% Last modified by fjsimons-at-alum.mit.edu, 06/13/2007

celnr=celnr(:);
celcum=[0 ; cumsum(nmr)];
rows=gamini(1:length(nmr),nmr');
rownr=rows(celnr)';
colnr=celnr-celcum(rownr);
if nargin>4
  c11=varargin{1};
  % Corresponding latitudes
  t=lat(rownr);
  b=lat(rownr+1);
  % Corresponding longitudes
  l=c11(1)+(colnr-1).*dlon(rownr);
  r=c11(1)+colnr.*dlon(rownr);
end
outs=[{ 't'} { 'b'} { 'l'} { 'r'}];
if nargout>2
  for index=1:(nargout-2)
    varargout{index}=eval(outs{index});
    end
end

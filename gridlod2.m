function [c11,cmn,lat,dlon,nmr,lej,dieptes,depth,disc,basdep,basfun,varargout]=...
    gridlod2(genx)
% [c11,cmn,lat,dlon,nmr,lej,dieptes,depth,disc,basdep,basfun,X,Y,XI,YI,refarea]=...
% GRIDLOD2(genx)
%
% Loads all the information pertaining to the three-dimensional grid
%
% INPUT:
%
% genx         string that identified the $GRIDS/directory with grid info
%
% OUTPUT:
%
% c11         [x,y]/[lon,lat]-coordinates of the upper left grid corner node [degrees]
% cmn         [x,y]/[lon,lat]-coordinates of the bottom right corner node [degrees]
% dlat        latitude interval, maintained throughout [degrees]
% dlon        longitude interval at the equator [degrees]
% nmr         the number of cells per row in the grid
% lej         the number of layers at depth <= until here contained in GRIDINFO
% dieptes     the vector with grid depth names as strings
% depth       the vector with grid depths as values <= until here loaded through GENERIC.M
% basdep      the central depth of the basis functions
% basfun      the actual basis function at depth <= until here contained in INAMAK3D
% X,Y         midpoint, i.e. pixel-centered, coordinates of the grid [degrees]
% XI,YI       midpoint coordinates of a regular interpolatory grid [degrees]
% refarea     reference area, size of one grid cell at the equator
%
% Last modified by fjsimons-at-alum.mit.edu, 08/20/2007

% Constructs the directory/filename containing the surface grid info
fn=fullfile(getenv('GRIDS'),genx,'gridinfo');

% Loads the basis functions identifying the third dimension of the grid
basfun=load(fullfile(getenv('GRIDS'),genx,'inamak3d'));

% Depth of the basis functions
basdep=basfun(:,1);
% The actual basis functions
basfun=basfun(:,2:end);

% Open the grid info file
fid=fopen(fn,'r');

% Loads grid corner points
fg=sscanf(fgetl(fid),'%f %f %f %f');
% c11 and cmn are NODE-CENTERED BOUNDARY POINTS unlike X, XI, Y and YI
c11=fg(1:2); c11=c11(:)';
cmn=fg(3:4); cmn=cmn(:)';

% Loads the remaining grid info
fg=sscanf(fgetl(fid),'%f %f %f');
% Latitude spacing
dlat=fg(1);
% Longitude spacing
dlon=fg(2);
% Number of depth layers
lej=fg(3);
fclose(fid);

% Make the equal-area grid based on this information
[lat,dlon,refarea,nmr]=authalic(c11,cmn,dlat,dlon);

% Run the 'generic' program that sets the depths and the discontinuities
addpath(fullfile(getenv('GRIDS'),genx))
generic
rmpath(fullfile(getenv('GRIDS'),genx))

if nargout>=13
  % Calculate midpoints of grid
  ranges=[zeros(size(nmr)) nmr-1]';
  ranges=matranges(ranges(:)');
  longes=gamini(dlon,nmr);
  lonran=ranges.*longes;
  latmp=[lat(1:end-1)+diff(lat)/2];
  lefts=gamini(c11(1)+dlon/2,nmr);
  % X, XI and Y, YI are PIXEL-CENTERED midpoints of (regular) grid
  X=lefts+lonran;
  Y=gamini(latmp,nmr);
  varargout{1}=X(:);
  varargout{2}=Y(:);
  if nargout==15
    XI=c11(1)+dlon(1)/2+[0:nmr(1)-1]*dlon(1);
    YI=latmp;
    varargout{3}=XI(:)';
    varargout{4}=YI(:);
  end
  if nargout==16
    varargout{5}=refarea;
  end
end

function varargout=plm2cube(lmcosi,lfin,alfa,bita,gama)
% v=PLM2CUBE(lmcosi,lfin,alfa,bita,gama)
%
% Transforms an array of real spherical harmonic coefficients into the
% faces of a cubed sphere.
%
% INPUT:
%
% lmcosi   The usual spherical harmonic array, as for PLM2XYZ
% lfin     Number of subdivisions [default: 6]
% alfa     First Euler angle of wholesale tilt of all tiles [defaulted]
% beta     Second Euler angle of wholesale tilt of all tiles [defaulted]
% gama     Third Euler angle of wholesale tilt of all tiles [defaulted] 
%
% OUTPUT:
%
% vp       The NxNx6 values to plot on the cubed sphere
%
% EXAMPLE:
%
% plm2cube('demo1','3D',1) % With topography and the continents in 3D
% plm2cube('demo1','2D',1) % With topography and the continents in 2D
% plm2cube('demo2','3D',1) % With Ritsema's model and the continents in 3D
% plm2cube('demo2','2D',1) % With Ritsema's model and the continents in 2D
%
% Last modified by fjsimons-at-alum.mit.edu, 2/23/2011

if ~isstr(lmcosi)
  % Default values
  defval('lfin',7)
  defval('alfa',[]);
  defval('bita',[]);
  defval('gama',[]);
  
  % Construct the cubed sphere coordinates
  [x,y,z]=cube2sphere(lfin,alfa,bita,gama);
  
  % If it takes a while, could expand the whole thing and interpolate to
  % the grid
  
  %h=waitbar(0);
  % Loop over the faces
  for in=1:6
    % waitbar(in/6,h,sprintf('Starting face %i/%i',in,6)) 
    % Transform these to longitude and latitude
    [phi,piminth,r]=cart2sph(x(:,:,in),y(:,:,in),z(:,:,in));
    lon=phi(:)*180/pi; lat=piminth(:)*180/pi;
    % Expand the spherical harmonics onto these grid points
    % Since I changed my mind on the parity, try both
    try
      v(:,:,in)=flipud(reshape(plm2xyz(lmcosi,lat,lon),...
			       [2^lfin+1 2^lfin+1]));
    catch
      v(:,:,in)=flipud(reshape(plm2xyz(lmcosi,lat,lon),...
			       [2^lfin 2^lfin]));
    end
  end
  % close(h)
  varns={v};
  varargout=varns(1:nargout);
elseif strcmp(lmcosi,'demo1')
  % This takes the place of 'tipe'
  defval('lfin','3D')
  % This takes the place of 'kind'
  defval('alfa',1)
  % Load Earth's topography
  lmcosi=fralmanac('GTM3AR','SHM');
  % Restrict to a manageable size
  L=48;
  lmcosi=lmcosi(1:addmup(L),:);
  % Perform the transform with standard inputs
  v=plm2cube(lmcosi,5);
  clf
  % Plot using "lfin" as an option
  [p,pg]=plotoncube(v,lfin,alfa);
  if strcmp(lfin,'3D')
    hold on; plotcont([],[],3)
    % Adjust the color scale to the data limits
    caxis(minmax(v))
  elseif strcmp(lfin,'2D')
    % Load and convert continents
    hold on
    [a,b,XYZ]=plotcont([],[],3); 
    delete(b)
    [phi,piminth,r]=cart2sph(XYZ(:,1),XYZ(:,2),XYZ(:,3));
    lonc=phi*180/pi; latc=piminth*180/pi;
    [xic,etac]=sphere2cube(lonc,latc);
    hold on
    [pc,pgc]=plotonchunk(xic,etac);
    hold off
  end
elseif strcmp(lmcosi,'demo2')
  % This takes the place of 'tipe'
  defval('lfin','3D')
  % This takes the place of 'kind'
  defval('alfa',1)
  % Load Ritsema's model at near 120 km depth
  %lmcosi=interpJRmodel(112.875,'xpcube');
  % Load Ritsema's model at near 410 km depth
  lmcosi=interpJRmodel(406.350,'xpcube');
  % Perform the transform with standard inputs
  %v=plm2cube(lmcosi);
  v=plm2cube(lmcosi,4);
  clf
  % Plot using "lfin" as an option
  [p,pg]=plotoncube(v,lfin,alfa);
  if strcmp(lfin,'3D')
    hold on; plotcont([],[],3)
    % Adjust the color scale to the data limits
    caxis(minmax(v))
  elseif strcmp(lfin,'2D')
    % Load and convert continents
    hold on
    [a,b,XYZ]=plotcont([],[],3); 
    delete(b)
    [phi,piminth,r]=cart2sph(XYZ(:,1),XYZ(:,2),XYZ(:,3));
    lonc=phi*180/pi; latc=piminth*180/pi;
    [xic,etac]=sphere2cube(lonc,latc);
    hold on
    [pc,pgc]=plotonchunk(xic,etac);
    hold off
  end
end

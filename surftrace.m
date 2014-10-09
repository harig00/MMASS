function surftrace(genx,fn,varargin)
% SURFTRACE(genx,fn,mkplot,rego)
%
% The main routine for the surface-wave azimuthal anisotropy inversion:
% surface-wave ray tracing and the constructions of sensitivity matrices
% for the inversion for isotropic heterogeneity and azimuthal anisotropy.
%
% INPUT:
%
% genx     string identifying $GRIDS directory with grid info for GRIDLOD2
% fn       string of the filename that contains event-station longitudes/latitudes
% mkplot   1 makes a plot as you go
% rego     0 for equal-area grids contained in the instructions in $GRIDS/genx
%          1 for regularly spaced grids instead
%
% OUTPUT:
%
% Output is to matrices C, L, S1 and S2 where 'S1' contains the path lengths
% and 'S2' the course of the tractory. 
%
% SEE ALSO: SURFDENS, SURFAZI
%
% Last modified by fjsimons-at-alum.mit.edu, 08/20/2007

% NOTE THIS MUST HAVE BEEN GENERATED FROM THE SAME
% DATA.NEWT FILE AS YOU'LL BE USING IN AMAK3d.

% Gets all the parameters of the three-dimensionsl grid
[c11,cmn,lat,dlon,nmr,lej,dieptes,depth,disc,basdep,basfun,X,Y,XI,YI,refarea]=...
    gridlod2(genx);

if nargin>=4 & varargin{2}==1
  disp('Making regular grid')
  nmr=repmat(max(nmr),length(nmr),1);
  dlon=repmat(min(dlon),length(dlon),1);
end

% Load stations and events
load(fn)
eval([ 'lonlat=' fn ';'])
% First of all, adhere to the convention 0<=lon<=360
lonlat(:,1)=lonlat(:,1)+(lonlat(:,1)<0)*360;
lonlat(:,3)=lonlat(:,3)+(lonlat(:,3)<0)*360;

% Check the boundedness of the grid
gridcheck(lonlat,c11,cmn)

[ms,ns]=size(lonlat);

% Open sensitivity matrix vector files
fid1=fopen('L.bin','w+');  % Cumulative nr of entries per row (event)
fid2=fopen('S1.bin','w+'); % Values of all entries (path lengths in km)
fid3=fopen('S2.bin','w+'); % Values of all entries (azimuths)
fid4=fopen('C.bin','w+');  % Cell (column) number of each entry

cumnum=0;
counter=0;
% Loop over paths
more off
for index=1:ms
  wb=waitbar(index/ms);

  lon1=lonlat(index,1);
  lat1=lonlat(index,2);
  lon2=lonlat(index,3);
  lat2=lonlat(index,4);

  % Calculate individual matrix elements
  [seqnr,gcd,crs,lonsec,latsec]=...
      gridsec(lat,dlon,[c11 cmn],[lon1 lat1],[lon2 lat2],nmr);
 
  % Make the plot as you go
  if nargin>=3 & varargin{1}==1
    clf
    fillauth(lat,dlon,nmr,c11,seqnr,gcd);    
    hold on
    plot(lonsec,latsec,'y-+')
    % plotcont([80 10],[220 -60])
    % axis([80 220 -60 10])
    axis image
    title(num2str(index),'FontSize',20)
    pause
  end
  
  fwrite(fid2,gcd,'float32');
  fwrite(fid3,crs,'float32');
  fwrite(fid4,seqnr,'int32');
  
  cumnum=[cumnum length(seqnr)];
  counter=counter+1;
  disp(sprintf('%4i %s',[counter 'events processed']))

end
more on
cumnum=cumsum(cumnum(2:end));

% Write last matrix
fwrite(fid1,cumnum,'int32');

% Close file
fclose('all');
close(wb)

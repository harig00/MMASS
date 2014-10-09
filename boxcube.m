function varargout=boxcube(lons,lats,deps,nxi,neta)
% BOXCUBE(lons,lats,deps,nxi,neta)
%
% Finds the coordinates, as a running index into the cubed-sphere, of a
% volume bounded two lines of longitude and two lines of latitude and two
% depths, or, alternatively, of a set of longitude, latitudes, and depths.
%
% INPUT:
%
% lons         A vector of longitudes (degrees), OR: two edge values
% lats         A same-size vector if latitudes (degrees), OR: two edge values
% deps         A vector of depths (m), OR: two edge values
% nxi,neta     Number of elements in the xi,eta direction
%
% OUTPUT:
%
% rinds         Vector with the running indices inside the cubed sphere
% xipi,etapi    A cell array with the face indices in xi and eta per chunk
% rindf         Cell array with the running indices on a chunk face
% rindc         Cell array with the running indices inside a chunk
% suming        The numbers of the faces that are nonempty in this respect
% dind          The numbers of the depth layers that are concerned
%
% EXAMPLE:
%
% boxcube('demo1') % Continental outlines on the faces of a cubed sphere
% boxcube('demo2') % An actual depth anomaly
% 
% See also: SPHERE2CUBE
%
% Last modified by fjsimons-at-alum.mit.edu, 01/27/2011

% Simple defaults
defval('lons',[100:110])

if ~isstr(lons)
  defval('lats',[20:30])
  defval('deps',[1])
  defval('nxi',128)
  defval('neta',128)

  % The spatial sampling interval on a face
  dxi=pi/2/(nxi-1);
  deta=pi/2/(neta-1);
  
  % Is it a box?
  if all([prod(size(lons)) prod(size(lats)) prod(size(deps))]==[2 2 2])
    lons=sort(lons); lats=sort(lats); deps=sort(deps);
    % Form all the pairs of points in a "square" patch, play it safe
    intv=180*min(dxi,deta)/pi/4;
    [lons,lats]=meshgrid(unique([lons(1):intv:lons(2) lons(2)]),...
			 unique([lats(1):intv:lats(2) lats(2)]));
    % The depth resolution is a hard constraint from the construction
    ddep=22.575;
    deps=unique([deps(1):ddep:deps(2) deps(2)]);
  end

  % Vectorize
  lons=lons(:)';
  lats=lats(:)';
  deps=deps(:)';

  % Check that they are pairs
  if size(lons)~=size(lats)
    error('Input must be pairs of longitudes/latitudes')
  end
  
  % This gives the xi,eta coordinates of the requested input as a cell
  % array per chunk, i.e. the face coordinates only
  [xip,etap]=sphere2cube(lons,lats);
  
  % Running index for depth, we do this first
  % Standard practice from PLOTONCUBE2
  % This here needs to be had from cubeunfolded.pdf
  load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS','radii_129'))
  % Now downsample this to 37 layers
  defval('depths',radii(129)-...
	 [radii([2 4 7:4:95 98 100 103 107 110 112 115:4:123 126])...
	  radii(127)+[radii(128)-radii(127)]/2 ...
	  radii(128)+[radii(129)-radii(128)]/2])
  
  % Nearest-neighbor to the depth, find depth index
  % Running index into a single chunk
  dind=unique(interp1(depths,1:length(depths),deps,'nearest'));

  % Figure out the index of xip and etap into the cubed-sphere array
  [xipi,etapi,rindf]=deal(cellnan(6,length(lons),1));
  suming=[];
  for chunk=1:6
    % Find the unique pairs of nearest neighbor
    xepu=unique(round(...
	[(xip{chunk}+pi/4)/dxi (etap{chunk}+pi/4)/deta])+1,'rows');

    % Here are the linear chunk indices that are hit
    xipi{chunk}=xepu(:,1);
    etapi{chunk}=xepu(:,2);

    % Running index on a face without depth
    rindf{chunk}=xipi{chunk}+(etapi{chunk}-1)*nxi;

    % Running indices over the depths in a single chunk
    rindc{chunk}=repmat(rindf{chunk},1,length(dind))+...
	repmat((dind-1)*nxi*neta,length(rindf{chunk}),1);
    
    % Running indices in the entire cubed sphere
    rinds{chunk}=rindc{chunk}+repmat((dind-1)*nxi*neta,...
				     size(rindc{chunk},1),1);
    
    % Figure out which ones are actually not empty
    if prod(size(xepu))~=0
      suming=[suming chunk];
    end

    % Linearize
    rindc{chunk}=rindc{chunk}(:);
    rinds{chunk}=rinds{chunk}(:);
  end
  
  % Convert to one running index for all by adding the chunks
  rinds=cat(1,rinds{:});
  
  % Output
  varns={rinds,xipi,etapi,rindf,rindc,suming,dind};
  varargout=varns(1:nargout);
  
elseif strcmp(lons,'demo1')
  tryme={'africa','australia','namerica','greenland','australia','eurasia'};
  ourpick=tryme{ceil(rand*length(tryme))};
  lonlats=eval(ourpick);
  N=2^round(3+rand*6);
  M=2^round(3+rand*6);
  display(sprintf('%s on a %ix%i cubed sphere',ourpick,N,M))
  [rinds,xipi,etapi,rindf,rindc,suming]=...
      boxcube(lonlats(:,1),lonlats(:,2),[],N,M);
  [X,fnX]=nancube(N,M,1);
  v=nan(N,M,6);
  for index=1:6
    X.(fnX{index})(rindf{index})=1;
    v(:,:,index)=X.(fnX{index});
  end
  [h,pg]=plotoncube(v,'2D',[],[],[],[],[],[],0);
  t(1)=title(sprintf(...
      'Outline of %s indexed to the voxels of a %ix%i cubed sphere',...
      upper(ourpick),M,N));
  xloc=-0.75;
  for chunk=1:length(suming)
    ril=length(rindf{suming(chunk)});
    t(1+chunk)=text(xloc,3.5-(chunk-1)*0.25,...
	      sprintf('chunk %i : %4i or %2i%s face indices',...
		      suming(chunk),ril,round(ril/M/N*100),'%'));
  end
  rilall=length(rinds);
  t(2+chunk)=text(xloc,3.5-chunk*0.25,...
		  sprintf('sphere  : %5i or %2i%s face indices',...
		      rilall,round(rilall/M/N/6*100),'%'));
elseif strcmp(lons,'demo2')
  % Standard 6 faces and 37 layers
  nxi=2^round(3+rand*4);
  neta=2^round(3+rand*4);
  ndeps=37;
  display(sprintf('Working on a %ix%i cubed sphere',nxi,neta))
  [rinds,xipi,etapi,rindf,rindc,suming,dind]=...
      boxcube([60 110],[-45 -5],[120 200],nxi,neta);
  % Assign these nonempties a certain value
  [X,fnX]=nancube(nxi,neta,ndeps);
  for index=1:6
    X.(fnX{index})(rindc{index})=1;
  end
  % Plot a single randomly chosen nonempty depth
  dshell=dind(ceil(rand*length(dind))); 
  [h,pg]=plotoncube2(X,dshell,'2D',[],[],[],0);
  xloc=-0.75;
  rilall=length(rinds); dimall=nxi*neta*ndeps*6;
  t=text(xloc,3.5-3*0.25,...
		  sprintf('sphere  : %5i / %7i or %3.1g%s face indices',...
		      rilall,dimall,rilall/dimall*100,'%'));
%  pause
%  plotoncube3(X,'2D',[],[],[],0);
%  keyboard
end

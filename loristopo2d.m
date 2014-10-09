function varargout=loristopo2d(wav,spp,N,J,L,precon,iface,colmap,tperc)
% [h,f,g,gg,t,tix,vw,bl]=loristopo2d(wav,spp,N,J,L,precon,iface,colmap,tperc)
%
% Plots a SINGLE-CHUNK wavelet-based multiscale rendition of TOPOGRAPHY
% on the cubed sphere using orthonormal wavelets. Called by LORIS3.
% 
% INPUT:
%
% wav        Wavelet basis ['D2' | D4' | 'D6']
% spp        0 You get a map of the truncated reconstruction
%            1 a plot of the truncated coefficients themselves
% N          The power of the dyadic subdivision [defaulted]
% J          Maximum scale in both directions [defaulted]
% L          Spherical harmonic degree limit [defaulted]
% precon     Two-array identifying preconditioners in ANGULARD{4,6}WT
% iface      The number of the special face that you are plotting
% colmap     Colormap string {'demmap' | 'kelicol' | 'sergeicol'}
% tperc      Truncation percentage (hard thresholding at percentil)
%
% OUTPUT:
%
% [h,f,g,gg,t] Handles to plot, grid, color bar, label, and title
% recerper     [1] Reconstruction error percentage
%              [2] numzero as a percentage of the total
%              [3] the actual truncation value that has been used 
%              [4] percentile absolute value that is truncated
% tix          The tick marks on the color bar and what percentiles they are
% vw           The original wavelet/scaling data
% bl           The thresholded reconstruced space-domain data
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/10/2012

% Set the defaults
defval('wav','D4')
defval('spp',1)
defval('N',7)
defval('J',N-4)
defval('L',ceil(2^(N+1)))
defval('precon',[1 1]*0)
defval('iface',3)
defval('colmap','kelicol')

% Color map saturation as percentiles
colperc=[5 95];
% For SPIE change it
colperc=[10 90];

% Truncation level as percentiles
defval('tperc',85);
% More percentiles for the tick marks on the color bar
ticperc=[50];

% Load or make the data, that is, the topography
fname=fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS',...
	       sprintf('loris2_%i_%i.mat',N,L));

if exist(fname,'file')==2
  load(fname)
else
  % Load Earth's topography
  lmcosi=rindeks(fralmanac('GTM3AR','SHM'),1:addmup(L));
  
  % Perform the transform with standard inputs
  v=plm2cube(lmcosi,N);

  % Save the data for the next time you run this
  save(fname,'v','N','L')
end

% This should be on pixel centered, non-overlapping, completely covering,
% registration with an even number of points. Right now we have this on
% grid node registration, with an odd number of points. Later we will
% change this, here we will fake this for the moment, to be quick
disp('Temporary faking of the pixel center registration')
v=v(1:2^N,1:2^N,:);

% Perform the thresholded reconstruction
[bl,vw,vwt,recerper,kilme]=angularthresh(v,tperc,J,iface,wav,precon);

% Now for the actual plot change the zeroes to NaN's
vwt(kilme)=NaN;

j=gca;
if spp==1
  % PLOT THE WAVELET AND SCALING COEFFICIENTS
  % Explicit and absolute color limits of the VALUE of the coeffs
  dax=prctile(vwt(:),colperc);
  h=imagefnan([1 1],[2^N 2^N],vwt,colmap,dax,[],[],0);
  hold on
  % Plot the wavelet grid
  f=fridplotw(N,J);
  hold off
  % Later redo the tick marks to more percentiles, e.g. the median
  tix=unique([dax prctile(vwt(:),ticperc)]);
  axis off
  strunk=sprintf('%s coefficients (m)',wav);
else
  % PLOT THE RECONSTRUCTED MAP VIEW
  % Explicit and absolute color limits of the VALUE of the recons
  dax=prctile(bl(:),colperc);
  h=imagefnan([1 1],[2^N 2^N],bl,colmap,dax,[],[],0);
  f=NaN;
  % Later redo the tick marks to more percentiles, e.g. the median
  tix=unique([dax prctile(bl(:),ticperc)]);
  noticks(gca)
  strunk=sprintf('%s reconstruction (m)',wav);
end
[g,gg]=addcb('hor',dax,dax,colmap,range(dax)/4);
% Note, as in CBARTICKS, that xlim and xtick aren't exactly the same
% as I am using IMAGEFDIR for the color bar - maybe time for a change
set(g,'Xtick',scale(tix,get(g,'xlim'))+...
      [1e-10 zeros(1,length(tix)-2) -1e-10],...
      'XtickL',round(tix))
axes(j)

stronk='rmse %3.1f%s | %i%sile | %i%s zero';
set(gg,'string',strunk)
t=title(sprintf(stronk,...
		   round(10*recerper(1))/10,'%',recerper(4),'%',...
		   round(recerper(2)),'%'));

% Make the output more informative
tix=[tix ; unique([colperc ticperc])];

% Optional output
varns={h,f,g,gg,t,recerper,tix,vw,bl};
varargout=varns(1:nargout);

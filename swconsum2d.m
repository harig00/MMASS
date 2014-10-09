function swconsum2d(par)
% SWCONSUM2D(par)
%
% FIGURE 6.2 of SIMONS & WANG
%
% Plots cumulative eigenvalue-weighted energy of the eigenfunction for
% various continental regions of concentration
%
% INPUT:
%
% par     1 Colorado plateaus [default]
%         2 Columbia plateau
%         3 Ozark plateaus
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2011

% Defaults
defval('par',1)
% What is the Shannon number?
defval('N',10)
% How many functions do you want?
defval('J',2*N)
% How many panels per row
defval('F',3)

% Load the USGS Tapestry
dirname=fullfile(getenv('IFILES'),'GEOLOGY','NORTHAMERICA');
load(fullfile(dirname,'tapestry.mat'))

% Switch the regions
names={'ColoradoPlateaus','ColumbiaPlateau','OzarkPlateaus'};

% Precomputed with defaults or not
fname=fullfile(getenv('IFILES'),'KERNELC2D',...
	       sprintf('SWREG-%s-%i-%i.mat',names{par},N,J));

if exist(fname,'file')==2
  load(fname)
  disp(sprintf('Loading %s',fname))
else
  % Find the curve defining the regions in GEOCENTRIC LONGITUDE and LATITUDE
  lola=[lon.(names{par})' lat.(names{par})'];
  % Figure out the center of mass for this projection
  [lonc,latc,A,lonm,latm]=rcenter(lola);
  % Check this with the SPHAREA computation which should be roughly the same
  difer(A-4*pi*spharea(lola),4)
  % Convert this to a 2D transverse-Mercator projection
  M=defaultm('tranmerc');
  M.origin=[latc lonc 0];
  % Use the WGS84 as the reference ellipsoid "datum" for the TM "projection"
  M.geoid=almanac('earth','wgs84','kilometers');
  % Fill in more stuff that is required
  M=defaultm(M);
  [X,Y]=mfwdtran(M,lola(:,2),lola(:,1));
  % Localize around the region with all the defaults
  [G,H,V,K,XYP,XY,A,c11,cmn]=localization2D([X Y],N,J);
  save(fname,'G','H','V','K','XYP','XY','A','c11','cmn')
end

% Perform the eigenvalue-weighted sum after initialization
SG=nan(size(G));
% Flipud
G=flipdim(G,1);
for in=1:J
  % Get the power spectrum of the tapers - normalized as in SWREGIONS2D
  SG(:,:,in)=V(in)*abs(fftshift(fft2(G(:,:,in)))).^2/prod(size(G(:,:,1)));
  G(:,:,in)=V(in)*G(:,:,in).^2;
end
% Check that Parseval stays the same
difer(sum(G(:))-sum(SG(:)))
% Perform progressive sums
elg(:,:,1)=sum(G(:,:,1:N/2),3);
elg(:,:,2)=sum(G(:,:,1:N),3);
elg(:,:,3)=sum(G(:,:,1:2*N),3);
kelg(:,:,1)=sum(SG(:,:,1:N/2),3);
kelg(:,:,2)=sum(SG(:,:,1:N),3);
kelg(:,:,3)=sum(SG(:,:,1:2*N),3);

% Derive a decent color scale for all plots
colscale=[-N/A N/A];
% The spectral color scale is not as simple - we work with spatially
% normalized eigenfunctions that are bandlimited so they all have power
% in the band and their complete sum should add up to A/4/pi^2 but since
% we start from the spatial functions we can never get there, and the
% truncated or eigenvalue-weighted sum is not good enough. Thus work in
% decibels as any normal individual would, anyways - see, e.g. PW 340-341.
kcolscale=[-40 0];

% Titles
tits={'N^{2D}/2' 'N^{2D}' '2N^{2D}'};

% Plot them
clf
[ah,ha,H]=krijetem(subnum(2,F));

defval('bw',1)

% In the space domain
for in=1:F
  axes(ah(in))
  % Check if the color saturation is appropriately guessed
  if bw==1
    imagefnan(c11,cmn,setnans(elg(:,:,in)),gray(21),colscale,grey(5))
  else
    imagefnan(c11,cmn,setnans(elg(:,:,in)),kelicol,colscale)
  end
  hold on
  p(in)=plot(XY(:,1),XY(:,2),'k');
  hold off
  xl(in)=xlabel('easting (km)');
  yl(in)=ylabel('northing (km)');
  tl(in)=title(sprintf('%s = 1 %s %s',...
		       '\alpha','\rightarrow',tits{in}));
  hold on
  plot([0 0],ylim,':')
  plot(xlim,[0 0],':')
  hold off
  switch par
   case 1
    set(gca,'xtick',[-500 0 500],'ytick',[-500 0 500])
   case 2
    set(gca,'xtick',[-500 0 500],'ytick',[-400 0 400])
  end
end

% Do you want a real wavelength scale or not?
defval('lambdas',0)

% In the spectral domain
for in=1:F
  axes(ah(in+F))
  % Make the frequency axis
  [fx,fy,fxnyq,fynyq,dx,dy]=...
      fftaxis(size(G(:,:,1)),size(SG),[range(XYP(:,2)) range(XYP(:,1))]);
  % The wavenumber axes corner points for after FFTSHIFT where you should
  % check that the symmetry point is properly at (0,0)
  kc11=[-max(fx) max(fy)]*2*pi;
  kcmn=[-min(fx) min(fy)]*2*pi;
  if ~lambdas
    kc11=kc11./[fxnyq fynyq]/2/pi;
    kcmn=kcmn./[fxnyq fynyq]/2/pi;
  end
  % Check if the color saturation is appropriate - don't switch DECIBEL and SETNANS
  if bw==1
      imagefnan(kc11,kcmn,decibel(setnans(kelg(:,:,in))),gray(21),kcolscale,grey(5))
  else
    imagefnan(kc11,kcmn,decibel(setnans(kelg(:,:,in))),kelicol,kcolscale)
  end
  hold on
  plot([0 0],ylim,':')
  plot(xlim,[0 0],':')
  % Plot the so-called bandwidth on there also
  % Note that the center is an off-centered pixel in the lower right
  % quadrant
  if lambdas
    xl(in+F)=xlabel('k_x (rad/km)');
    yl(in+F)=ylabel('k_y (rad/km)');
    fpp(in)=circ(K);
    axis image
  else
    xl(in+F)=xlabel('k_x/k_x^{Nyq}');
    yl(in+F)=ylabel('k_y/k_y^{Nyq}');
    % Note this assumes dx==dy to normalize K
    fpp(in)=circ(K/pi*dx);
    % Only show some tiny fraction of this range
    axis([-1 1 -1 1]/10)
    set(gca,'xtick',[-0.1 0 0.1],'ytick',[-0.1 0 0.1])
  end
  hold off
end

% Cosmetics
nolabels(ha(3:end),2)
longticks(H)
delete(yl(2:F))
delete(yl(F+2:2*F))
delete(ah(2*F+1:end))
fig2print(gcf,'portrait')
shrink(ah,1.2,1.2)
serre(H,[],'across')
movev(ah(F+1:end),.03)
movev(xl(F+1:end),-0.0075)
if par==2
  movev(ah(1:F),-.035)
end

% Need color scales
h=axes;
if bw==1
  [cb(1),xcb(1)]=addcb('hor',colscale,colscale,gray(21));
else
  [cb(1),xcb(1)]=addcb('hor',colscale,colscale,'kelicol');
end
delete(h)
shrink(cb(1),2,2); movev(cb(1),0.35)
set(xcb(1),'String','cumulative energy')
set(cb(1),'xlim',[0 colscale(2)])
set(cb(1),'xtick',[0 colscale(2)],'xtickl',{'0' 'N/A'})
movev(xcb(1),3)

h=axes;
if bw==1
  [cb(2),xcb(2)]=addcb('hor',kcolscale,kcolscale,gray(21));
else
  [cb(2),xcb(2)]=addcb('hor',kcolscale,kcolscale,'kelicol');
end
delete(h)
shrink(cb(2),2,2); movev(cb(2),-0.09)
set(xcb(2),'String','cumulative energy (dB)')
set(cb(2),'xlim',[-20 0])
set(cb(2),'xtick',[-20:5:0],'xtickl',[-20:5:0])

% Set both of these in the middle
hh=getpos(ah(2));
cc=hh(1)+hh(3)/2;

layout(cb',cc,'m','x')

figdisp([],par,[],0)


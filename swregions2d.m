function swregions2d(par)
% SWREGIONS2D(par)
%
% FIGURE 6.1 of SIMONS & WANG
%
% Plots spatial and spectral renditions of the basis functions
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
% How many functions are you plotting?
defval('F',5)

% Load the USGS Tapestry in GEOCENTRIC coordinates
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
  [lonc,latc,A2,lonm,latm]=rcenter(lola);
  % Check this with the SPHAREA computation which should be roughly the same
  difer(A2-4*pi*spharea(lola),4)
  % Convert this to a 2D transverse-Mercator projection
  M=defaultm('tranmerc');
  M.origin=[latc lonc 0];
  % Use the WGS84 as the reference ellipsoid "datum" for the TM "projection"
  M.geoid=almanac('earth','wgs84','kilometers');
  % Fill in more stuff that is required
  M=defaultm(M);
  % The output of this is now in km and thus the area in square kilometers
  [X,Y]=mfwdtran(M,lola(:,2),lola(:,1));
  % Localize around the region with all the defaults
  [G,H,V,K,XYP,XY,A,c11,cmn]=localization2D([X Y],N,J);
  % Could check the area again
  disp(sprintf('Area determined on sphere   %8.0f km^2',...
	       A2*(fralmanac('Radius')/1000)^2))
  disp(sprintf('Area on projected ellipsoid %8.0f km^2',A))
  save(fname,'G','H','V','K','XYP','XY','A','c11','cmn')
end

% Derive a decent color scale for all plots
colscale=2*[-sqrt(1/A) sqrt(1/A)];

% This turns out to be a decent guess but it is 
kcolscale=8*[-1/pi/K.^2 1/pi/K.^2]/(4*pi^2)/prod(size(G(:,:,1)));

% Plot them
clf
[ah,ha,H]=krijetem(subnum(3,F));

defval('bw',1)

% In the space domain
for in=1:F
  axes(ah(in))
  % Check if the color saturation is appropriately guessed from the
  % orthogonality condition and the eigenvalue-weighted sum
  if bw==1
    imagefnan(c11,cmn,setnans(flipud(G(:,:,in))),...
              gray(21),colscale,grey(5))
  else
    imagefnan(c11,cmn,setnans(flipud(G(:,:,in))),...
              kelicol,colscale)
  end
  hold on
  p(in)=plot(XY(:,1),XY(:,2),'k');
  hold off
  tlb(in)=title(sprintf('%s_{%i} = %9.6f','\lambda',in,V(in)));
  xl(in)=xlabel('easting (km)');
  yl(in)=ylabel('northing (km)');
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
  % Get the data
  theG=flipud(G(:,:,in));
  % Get the power spectrum of the tapers - normalize
  SG=abs(fftshift(fft2(theG))).^2/prod(size(theG));
  % Check Parseval's relation - see also PARSEVAL
  difer(sum(theG(:).^2)-sum(SG(:)))
  % Make the frequency axis
  [fx,fy,fxnyq,fynyq,dx,dy]=...
      fftaxis(size(theG),size(SG),[range(XYP(:,2)) range(XYP(:,1))]);
  % The wavenumber axes corner points for after FFTSHIFT where you should
  % check that the symmetry point is properly at (0,0)
  kc11=[-max(fx) max(fy)]*2*pi;
  kcmn=[-min(fx) min(fy)]*2*pi;
  if ~lambdas
    kc11=kc11./[fxnyq fynyq]/2/pi;
    kcmn=kcmn./[fxnyq fynyq]/2/pi;
  end
  % Check if the color saturation is appropriately guessed from the
  % orthogonality condition and the eigenvalue-weighted sum and the lack
  % of normalization of the DFTMTX and the area in space vs the area in
  % the wavenumber space
  if bw==1
    imagefnan(kc11,kcmn,setnans(SG),gray(21),kcolscale,grey(5))
  else
    imagefnan(kc11,kcmn,setnans(SG),kelicol,kcolscale)
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
nolabels(ha(3:3*F),2)
longticks(H)

if lambdas
  % Wavelength axis needs to come next to last
  lambda=[-80 -20 -10 10 20 80];
  for in=1:F
    [ax(in),xlx(in),ylx(in)]=xtraxis(ah(in+F),...
		 2*pi./lambda,abs(lambda),'wavelength (km)',...
		 2*pi./lambda,abs(lambda),'wavelength (km)');
    rottick(ax(in),'t')
    movev(xlx(in),.25)
  end
  nolabels(ax(1:F-1),2)
  delete(ylx(1:F-1))
  longticks(ax)
  %grat=size(theG,1)/size(theG,2);
  set(ax,'PlotBoxA',get(ah(F+1),'PlotBoxA'))
  shrink(ax,1.025,1.025)
  moveh(ax,0.001)
else
  switch par
   case 1
    serre(H',2/3,'down')
    movev(xl(F+1:end),-0.01*0.8/2)
    moveh(yl(F+1:end),0.0025*0.8*15)
   case 2
    serre(H',1,'down')
    movev(xl(F+1:end),-0.01)
    moveh(yl(F+1:end),0.0025)
   case 3
    serre(H',1,'down')
    movev(xl(F+1:end),-0.01)
    moveh(yl(F+1:end),0.0025)
  end
  shrink(ah,1/1.5,1/1.5)
  movev(ah(1:F),0.01)
end

set([tlb xl(~~xl) yl(~~yl) ah(~~ah)],'FontS',12)
delete(yl(2:F))
delete(yl(F+2:2*F))
delete(ah(2*F+1:end))
fig2print(gcf,'landscape')
figdisp([],par,[],0)

keyboard

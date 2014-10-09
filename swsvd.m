function swsvd(par,tol)
% SWSVD(par,tol)
%
% FIGURE 6.3 of SIMONS & WANG
%
% Illustration of arbitrary Cartesian Slepian functions with fan-like
% anisotropic spectral sensitivity.
%
% INPUT:
%
% par     1 Colorado plateaus [default]
%         2 Columbia plateau
%         3 Ozark plateaus
% tol     abs(log10(tolerance)) for EIGS in SVDSLEP3 
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2011

% How many functions do you want?
defval('J',20)
% Defaults
defval('par',1)
defval('tol',17)
% How many panels are you plotting on one row?
defval('F',5)
% How big does the computation domain get?
defval('ngro',3)

% Load the USGS Tapestry
dirname=fullfile(getenv('IFILES'),'GEOLOGY','NORTHAMERICA');
load(fullfile(dirname,'tapestry.mat'))

% Switch the regions
names={'ColoradoPlateaus','ColumbiaPlateau','OzarkPlateaus'};

% Precomputed with defaults or not
fname=fullfile(getenv('IFILES'),'KERNELC2D',...
	       sprintf('SWSVD-%s-%i-%i.mat',names{par},tol,ngro));

if exist(fname,'file')==2 & 1==3
  load(fname)
  disp(sprintf('Loading %s',fname))
else
  % This here is from SWREGIONS2D but most of it gets thrown out
  N=10;
  load(fullfile(getenv('IFILES'),'KERNELC2D',...
	       sprintf('SWREG-%s-%i-%i.mat',names{par},N,J)));
  % We just want the area A and the outline XY
  clear G H V K XYP c11 cmn
  
  % Define a localization interval in spectral space
  defval('R',300);
  % A half triangle in spectral space
  KXY=sqrt(pi*(R/5)^2)/2*[0  1/2 -1/2 0;...
    		          0  1  1 0]';
  % Clockwise rotation i the Fourier domain!
  rr=rotz(pi/5);
  KXY=[rr(1:2,1:2)*KXY']';

  % Localize around the region with all the defaults
  [E1,V1,c11cmnR1,c11cmnK1,SE1,KXY1]=svdslep3(XY,KXY,J,tol,ngro);
  
  % A half triangle in spectral space
  KXY=sqrt(pi*(R/5)^2)/2*[0  1/2 -1/2 0;...
    		          0  1  1 0]';

  % Clockwise rotation i the Fourier domain!
  rr=rotz(-pi/5);
  KXY=[rr(1:2,1:2)*KXY']';

  % Localize around the region with all the defaults
  [E2,V2,c11cmnR2,c11cmnK2,SE2,KXY2]=svdslep3(XY,KXY,J,tol,ngro);

  % Save them all now 
  save(fname,'E1','V1','c11cmnR1','c11cmnK1','SE1','KXY1','XY',...
       'E2','V2','c11cmnR2','c11cmnK2','SE2','KXY2','A','ngro')
end

% Now make the figure
clf
[ah,ha,H]=krijetem(subnum(3,F));
axshr=1.5;
axshrK=4.75;

% Derive a decent color scale for space plots - compare SWREGIONS2d
colscale=10*[-sqrt(1/A) sqrt(1/A)];
% Derive a decent color scale for spectral plots
colscaleK=[-200 0];

% First the eigenfunctions
for ind=1:F-1
  axes(ah(ind))
  [tlb(ind),xl(ind),yl(ind)]=...
      plotspace(ind,E1,c11cmnR1,V1,axshr,XY,colscale);

  axes(ah(ind+F))
  [tlb(ind+F),xl(ind+F),yl(ind+F)]=...
      plotspace(ind,E2,c11cmnR2,V2,axshr,XY,colscale);
end

% Then the cumulative eigenvalue-weighted spectral energy
xte=c11cmnK2(3)-c11cmnK2(1);
yte=c11cmnK2(2)-c11cmnK2(4);


axes(ah(F))
[p(1),tlb(F),xl(F),yl(F)]=...
    plotspec(SE1,V1,c11cmnK1,KXY1,axshrK,colscaleK);

axes(ah(2*F))
[p(2),tlb(2*F),xl(2*F),yl(2*F)]=...
    plotspec(SE2,V2,c11cmnK2,KXY2,axshrK,colscaleK);

upit=22.5;

% Cosmetics
switch par
 case 1
  serre(H',2/3,'down')
  set(ha(1:3*F-4),'xtick',[-500 0 500],'ytick',[-500 0 500])
 case 2
  serre(H',1,'down')
  set(ha(1:3*F-4),'xtick',[-500 0 500],'ytick',[-400 0 400])
 case 3
  serre(H',1,'down')
end
set(ah([F 2*F]),'xtick',[-0.1 0 0.1]*xte,'xtickl',[-0.1 0 0.1],...
		'ytick',[-0.1 0 0.1]*yte,'ytickl',[-0.1 0 0.1])
sro=1.5;
shrink(ah,1/sro,1/sro)
movev(ah(1:F),0.03)
movev(tlb([F 2*F]),upit)
moveh(ha(1:3*F-4),-0.02)
nolabels(ah(1:F),1)
nolabels(ha(3:3*F-4),2)
set(ah([F 2*F]),'yaxisl','r')
longticks(H)
set([tlb(~~tlb) xl(~~xl) yl(~~yl) ah(~~ah)],'FontS',12)
delete(xl(1:F))
delete(yl(2:F-1))
delete(yl(F+2:2*F-1))
delete(ah(2*F+1:end))

axes(ah(2*F))
pos=[0.8005    0.26    0.1249    0.0275];
defval('bw',1)
if bw==1
  [cb,xcb]=addcb(pos,colscaleK,colscaleK,gray(21),25);
else
  [cb,xcb]=addcb(pos,colscaleK,colscaleK,'kelicol',25);
end

set(cb,'xlim',[colscaleK(1)/2 colscaleK(2)+1])
longticks(cb)

shrink(cb,1,1.5)
set(xcb,'string','cumulative energy (dB)')

fig2print(gcf,'landscape')
figdisp([],par,[],0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x,y]=plotspace(ind,E,c11cmn,V,axshr,XY,colscale)
toplot=v2s(E(:,ind));
defval('bw',1)
if bw==1
  imagefnan(c11cmn(1:2),c11cmn(3:4),toplot,gray(21),colscale,grey(5))
else
  imagefnan(c11cmn(1:2),c11cmn(3:4),toplot,[],colscale)
end

t=title(sprintf('%s_{%i} = %9.6f','\lambda',ind,V(ind)));
x=xlabel('easting (km)');
y=ylabel('northing (km)');
axis([c11cmn(1) c11cmn(3) c11cmn(4) c11cmn(2)]/axshr)
hold on
plot(XY(:,1),XY(:,2),'k'); 
plot([0 0],ylim,':')
plot(xlim,[0 0],':')
hold off
movev(x,-40)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,t,x,y]=plotspec(SE,V,c11cmn,KXY,axsh,colscale)
SEE=sum(repmat(V(:)',size(SE,1),1).*SE.^2,2);
psdens=fftshift(decibel(v2s(SEE)));
psdens(psdens<colscale(1))=NaN;
defval('bw',1)
if bw==1
  imagefnan(c11cmn(1:2),c11cmn(3:4),psdens,gray(21),colscale); 
else
  imagefnan(c11cmn(1:2),c11cmn(3:4),psdens,[],colscale); 
end
axis image
axis([c11cmn(1) c11cmn(3) c11cmn(4) c11cmn(2)]/axsh)
hold on
p=plot(KXY(:,1),KXY(:,2),'k');
plot([0 0],ylim,':')
plot(xlim,[0 0],':')
hold off
t=title(sprintf('%s = 1 %s %i',...
		     '\alpha','\rightarrow',length(V)));
x=xlabel('k_x/k_x^{Nyq}');
y=ylabel('k_y/k_y^{Nyq}');
movev(x,-5)

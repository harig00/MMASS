function dscmb
% Makes a plot of the CMB spectrum, its WS error, and how it is improved
% using eigenvalue-weighted multitaper analysis

% PARAMETERS:
% $\Omega_{\mathrm b}=0.046,
%\Omega_{\mathrm c}=0.224, \Omega_{\Lambda}=0.730$ and $H_0=72$
% $\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{Mpc}^{-1}$

% This is the file name where the spectrum is kept
fname='cmb_68217229.fcl';

% These values are alreay prewhitened, and in microKelvin^2
cmb=load(fullfile(getenv('IFILES'),'CMB',fname));

% Trim this a bit
xmax=900;
cmb=cmb(cmb(:,1)<=xmax,:);

% The degrees at which this spectrum is defined
l=cmb(:,1);
% CMB temperature in microKelvin
Tcmb=2.725e6;
% Whitened spectrum which contains the l(l+1)/2/pi factor
S=cmb(:,2).*Tcmb^2;

% Noise parameters of the detector
dom=4e-6;
% In microkelvin per pixel
sig=100;
% fwhm in arcmin now in radians now
th=20*(pi/180/60);
% Noise spectrum
N=dom*sig^2.*exp(l.^2*th^2/8/log(2));
% Whitened spectrum using the same whitening as th spectrum
N=l.*(l+1)/2/pi.*N;

% Theoretical standard error for the whole-sphere no-taper result
tstd=sqrt(2./(2*l+1)).*(S+N);

% DOUBLE CAP
sord=2;
% SIZE OF THE CUT
TH=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
[ah,ha]=krijetem(subnum(3,2));

xlabs=sprintf('degree l');
xlobs=sprintf('angular scale');
ylabs=sprintf('l(l+1)S_l/2%s (%sK^2)','\pi','\mu');

axes(ah(1))

EL=[10 20 30 40 50 60];

% Make the small panel plots
for inx=1:length(ah)
  [psl(inx),perr(inx),psb(inx),xl(inx),yl(inx)]=...
      doitsmoothit(EL(inx),cmb,S,l,tstd,ah(inx),xlabs,ylabs,sord,TH,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cosmetics
set(psl,'Color','k','LineW',1)
set(psb,'LineS','none','Marker','o','MarkerS',3,...
	'MarkerF','w','MarkerE','k')
xtix=[2 100:200:xmax];
set(ah,'xlim',[0 xmax],'ylim',[0 6000],'xtick',[2 100:200:xmax])
longticks(ah)
set(ah(~~ah),'box','on')
% Do this before the serre operation
delete(yl([2 4 6]))
delete(xl(1:4))
nolabels(ah(1:4),1)
nolabels(ha(4:6),2)

serre(ha(1:3),2/5,'down')
serre(ha(4:6),2/5,'down')

serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ah(5:6),1/2,'across')

for inx=1:length(ah)
  % Put the identifying labels on
  [bh,th]=boxtex('um',ah(inx),...
		 sprintf('L = %i',EL(inx)),12,1,0.75,1.15);
  % Put a caplogo on here
  cl(inx)=caplogo(ah(inx),sord+(sord==2),...
			     'ul',1/30,1/30);

end
movev(cl,-0.035)
moveh(cl,-0.002)
% Put the TH of the left-over cap on here
for inx=1:length(ah)
  axes(ah(inx))
  tx(inx)=text(10,4800,sprintf('%s = %i%s','\Theta',...
	       90*(sord==2)+(-1)^(sord+1)*TH,str2mat(176)));
end
set(tx,'FontS',8)

% Put an extra axis on top here using the Jeans relation
degs=[0.2 0.3 0.5 1 90];
for inx=1:length(degs)
  xtox(inx)=max(roots([1 1 -(180/degs(inx))^2]));
end
[ax(1),xxl(1)]=xtraxis(ah(1),xtox,degs,xlobs);
[ax(2),xxl(2)]=xtraxis(ah(2),xtox,degs,xlobs);
longticks(ax)
deggies(ax,1)
% On the others put an "empty" xtraxis to properly close the plots
for inx=3:6
  ax(inx)=xtraxis(ah(inx));
  set(ax(inx),'box','on')
end

% Put the caplogo back on top where it belongs
for inx=1:length(ah)
  axes(cl(inx))
end
  
set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')
figdisp([],[],[],0)
!degs /home/fjsimons/EPS/dscmb.eps

disp('Rerun to get dimensions exactly right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SB,B]=smoothit(Lmax,L,S)
% Calculate coupling matrices
% The maximum of the WIGNER0JC database
dbmax=500;
M=mcouplings(L,Lmax,1,dbmax-L);
% Take out rows so the bins become non-overlapping
% These are the bins, add one to get the bins indices
B=L/2:2*L+1:Lmax;
% Take out the degree 0 and 1 part from the true spectrum
M=M(B+1,3:end);
% Perform the deliberate biasing
SB=M*S;
hold on
% Fix the last element up for cosmetic reasons
SB(end)=SB(end)/sum(M(end,:),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psl,perr,psb,xl,yl]=...
    doitsmoothit(L,cmb,S,l,tstd,ah,xlabs,ylabs,sord,TH,xver)
axes(ah)
% Calculate the biased measurements
[SB,B]=smoothit(cmb(end,1),L,S);
% Plot the true values and the shading
% Make end bits a bit nicer by introducing fake data
% Note that this is NOT a condifence interval; we don't know the
% distribution really
[psl,e68,e95]=errorfill...
    ([l(1)-10 ; l(:) ; l(end)+10],[S(1) ; S(:) ; S(end)],...
     [tstd(1) ; tstd(:) ; tstd(end)]);
delete(e95)
set(psl,'Color','k','LineW',1)
set(e68,'FaceC',grey(9))

% Plot the truth without the shaded bars
% psl=plot(l,S);
hold on
% No, hold on, calculate based on a cut
[mt2ws,ml,mt2wsinf]=mvarratios(L,TH,sord,B,xver);
if xver==1
  % Calculate the error bars based on full coverage
  mt2wsF1=mvarratios(L,180,1,B,xver);
  mt2wsF2=mvarratios(L,0,2,B,xver);
  difer(mt2wsF1-mt2wsF2,[],[],'All consistent')
  % If you used a ratio, it would come out
  % If you didn't, it shows you how it would have come out
  % We only show the case where it has been used, i.e. where it comes out
  % constant 
  rus=mt2ws./mt2wsF1;
  if ~sum(diff(rus))
    switch sord
     case 1
      Ato4pi=(1-cos(TH*pi/180))/2;
     case 2
      Ato4pi=(1-sin(TH*pi/180));
    end
    disp(sprintf('Exponent that appears to be used is %4.2f',...
		 log10(rus(1))/log10(Ato4pi)))
    disp(sprintf('Actual ratio appearing to be used is %4.2f',rus(1)))
  end
end
% Sure, at these bin widths, you may hit the maximum variance ratio first
% and then come down - you'd have to compute it at B=l to see the whole
% pattern 
disp(sprintf('Error ratio at large l is %5.3f',sqrt(mt2wsinf)))
% Plot these as error bars
perr=errorbar(B,SB,sqrt(mt2ws(:)).*tstd(B+1));
set(getkids(perr,2),'Color','k')
delete(getkids(perr,1))
hold on
% Plot the biased measurements
psb(1)=plot(B,SB);
yl=ylabel(ylabs);
xl=xlabel(xlabs);

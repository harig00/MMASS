function varargout=admittb(infl,onfl,land)
% ADMITTB(infl,onfl,land)
%
% Expected admittance for top loading and for bottom loading where the
% level of bottom loading zl is equal to compensating interface zm.
%
% Forsyth (1985), Figure 1.
%
% OUTPUT:
%
% Lots of handle variables
%
% See also ADMITF2, ADMITR, COHERETE, COHEREF2
%
% Last modified by fjsimons-at-alum.mit.edu, 03/07/2012

defval('infl',1)
defval('onfl',1)

% Earth model
rc=2670;
rm=3300;
drho=rm-rc;
g=fralmanac('GravAcc');
G=fralmanac('GravCst');

% Elastic parameters 
young=1.4e11; % [Nm-2]
poisson=1/4;

% Compensation depth % [m]
zm=35000;

% Elastic thicknesses [m]
Te=[0 5 10 20:20:100]*1000;

% Flexural rigiditiy [Pa*m^3 or N*m]
D=Te.^3*young/12/(1-poisson^2);
% Wavelengths
lambda=linspace(1,5000,2000)*1000;

[LL,DD]=meshgrid(lambda,D);
KK4=(2*pi./LL).^4;

% Forsyth Eqs. (3) and (6)
xai=1+DD.*KK4/drho/g;
phi=1+DD.*KK4/rc/g;
upco=exp(-(2*pi./LL)*zm);

% Top- and bottom-only admittances
% Olhede and Simons Eqs. (61) and (62)
% Could write this using ADMITTANCE
QT=-2*pi*G*rc*upco./xai;
QB=-2*pi*G*rc*upco.*phi;

% Convert to mgal/m
QT=QT/1e-5;
QB=QB/1e-5;

% Cross check with ADMITTANCE which itself has been cross-checked with FORSYTH
% for uncorrelated loading only
for in=1:length(D)
  difer(QT(in,:)-admittance(Te(in),0,0,zm,rc,drho,[],2*pi./lambda,0),...
	[],[],NaN);
  difer(QB(in,:)-admittance(Te(in),Inf,0,zm,rc,drho,[],2*pi./lambda,0),...
	[],[],NaN);
end

% Do the plotting proper
figure(gcf)
for index=1:length(D)
  pt(index)=semilogx(lambda/1000,QT(index,:),'k-');  
  hold on
end
for index=1:length(D)
  pb(index)=semilogx(lambda/1000,QB(index,:),'k--');  
  hold on
end
hold off

% Cosmetics 
set(gca,'xdir','rev','ydir','rev')
yls=[-0.2 0.01];
ylim(yls)
xlim([10 5000])

grid on
longticks

% xl=xlabel('Wavelength (km)');
% yl=ylabel('Bouguer per Topography: Admittance (mgal/m)');
xl=xlabel('wavelength (km)');
yl=ylabel('admittance Q_1, Q_2 (mgal/m)');

set(pt,'linew',2)
set(pb,'linew',2)

% Legend boxes
post=[0.6700    0.0900  % 0
      0.5900    0.0900  % 5 
      0.5072    0.1400  % 10
      0.4200    0.1800  % 20
      0.3400    0.1800  % 40
      0.2800    0.2300  % 60
      0.2400    0.2800  % 80
      0.2000    0.3300];% 100

posb=[0.2500    0.7200 % 100
      0.2800    0.6700 % 80 
      0.3100    0.6200 % 60
      0.3500    0.5700 % 40
      0.4800    0.5700 % 20
      0.6500    0.5700 % 10
      0.6657    0.1700 % 5
      0.040     0.0465];

hold on
b=gca;
d=axes;
tbb=0.015*infl;
tbx=0.03*onfl;
fb=fillbox(...
    [post(:,1)-tbx post(:,1)+tbx post(:,2)+tbb post(:,2)-tbb],'w');
hold on
fbb=fillbox(...
    [posb(1:end-1,1)-tbx posb(1:end-1,1)+tbx posb(1:end-1,2)+tbb posb(1:end-1,2)-tbb],'w');
fbb=[fbb ; fillbox(...
    [posb(end,1)-tbx/onfl posb(end,1)+tbx/onfl posb(end,2)+tbb*1.25 posb(end,2)-tbb*1.25],'w')];
axes(d)
axis([0 1 0 1])
set(d,'color','none')
nolabels(d)
noticks(d)

for index=1:length(Te)
  t(index)=text(post(index,1),post(index,2),...
		num2str(Te(index)/1000));
end
for index=1:length(Te)-1
  tb(index)=text(posb(index,1),posb(index,2),...
		 num2str(Te(length(Te)-index+1)/1000));
end
tb(index+1)=text(posb(index+1,1),posb(index+1,2),'Te');

set([t tb],'horizontala','center')

set(fb,'facecolor',[0.75 0.75 0.75])

set(d,'position',get(b,'position'))

set(b,'XTickLabel',[10 100 1000])
set([gca xl yl],'fonts',12)

defval('land',1)
if land==1
  fig2print(gcf,'landscape')
  figdisp
end

varns={d,xl,yl,pt,pb,t,tb};
varargout=varns(1:nargout);


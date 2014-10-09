function varargout=admitf2(infl,onfl,land)
% ADMITF2(infl,onfl,land)
%
% Expected admittance for top and bottom loading combined. Second
% compensating interface is same as level of loading. In function of
% bottom-to-top loading fraction f2.
%
% Forsyth (1985), Figure 5.
%
% OUTPUT:
%
% Lots of handle variables
%
% See also: ADMITR, ADMITB, COHERETE, COHEREF2
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

% Elastic thickness [m]
Te=40*1000;
% Loading ratios
f2=[0 0.5 1 2 5];

% Flexural rigiditiy [Pa*m^3 or N*m]
D=Te^3*young/12/(1-poisson^2);
% Wavelengths [m]
lambda=linspace(1,5000,2000)*1000;

[LL,FF2]=meshgrid(lambda,f2);
KK4=(2*pi./LL).^4;

% Forsyth Eqs. (3) and (6)
xai=1+D.*KK4/drho/g;
phi=1+D.*KK4/rc/g;
upco=exp(-(2*pi./LL)*zm);

% Forsyth Eq. (12)
BtoT=sqrt(FF2)*rc./xai/drho;

% Combined admittance for independent surface
% and subsurface loading
% Olhede and Simons Eq. (60), same thing
% Could write this using ADMITTANCE
QTB=-2*pi*G*rc*upco.*(phi.*BtoT.^2+1./xai)./(BtoT.^2+1);
% Rewrite another way just for fun
fax=sqrt(FF2)*rc./drho;
QTBR=-2*pi*G*rc*upco.*(xai+fax.^2.*phi)./(xai.^2+fax.^2);
difer(QTB-QTBR,[],[],NaN)

% Convert to mgal/m
QTB=QTB/1e-5;

% Cross check with ADMITTANCE which itself has been cross-checked with FORSYTH
% for uncorrelated loading only
for in=1:length(f2)
  difer(QTB(in,:)-admittance(Te,f2(in),0,zm,rc,drho,[],2*pi./lambda,0),...
	[],[],NaN);
end

% Do the plotting proper
figure(gcf)
for index=1:length(f2)
  pt(index)=semilogx(lambda/1000,QTB(index,:),'k-');  
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

%xl=xlabel('Wavelength (km)');
%yl=ylabel('Bouguer per Topography: Admittance (mgal/m)');
xl=xlabel('wavelength (km)');
yl=ylabel('admittance Q_f (mgal/m)');

set(pt,'linew',2)

% Legend boxes
pos=[0.3500    0.2300
     0.3700    0.2800
     0.3900    0.3300
     0.4100    0.3800
     0.4300    0.4300
     0.0375    0.0465];

hold on
b=gca;
d=axes;
tbb=0.015*infl;
tbx=0.03*onfl;
fb=fillbox(...
    [pos(1:end-1,1)-tbx pos(1:end-1,1)+tbx pos(1:end-1,2)+tbb pos(1:end-1,2)-tbb],'w');
hold on
fb=[fb ; fillbox(...
    [pos(end,1)-tbx/onfl pos(end,1)+tbx/onfl pos(end,2)+tbb*1.25 pos(end,2)-tbb*1.25],'w')];
axes(d)
axis([0 1 0 1])
set(d,'color','none')
nolabels(d)
noticks(d)

for index=1:length(f2)
  t(index)=text(pos(index,1),pos(index,2),num2str(f2(index)));
end
t(index+1)=text(pos(index+1,1),pos(index+1,2),'f^2');

set(t,'horizontala','center')

set(d,'position',get(b,'position'))

set(b,'XTickLabel',[10 100 1000])
set([gca xl yl],'fonts',12)

defval('land',0)
if land==1
  fig2print(gcf,'landscape')
  figdisp
end

varns={d,xl,yl,pt,t};
varargout=varns(1:nargout);

function varargout=admitr(infl,onfl,land)
% ADMITR(infl,onfl,land)
%
% Expected admittance for correlated top and bottom loading
% combined. Second compensating interface is same as level of loading. In
% function of bottom-to-top initial-loading ratio r.
%
% See Crosby (2007)
%
% OUTPUT:
%
% Lots of handle variables
%
% See also: ADMITF2, ADMITB, COHERETE, COHEREF2
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
% Loading ratio
f2=1;
% Correlation coefficients (r and -r and -1 will be considered)
r=[0 0.25 0.5 0.75];

% Flexural rigidity [Pa*m^3 or N*m]
D=Te^3*young/12/(1-poisson^2);
% Wavelengths [m]
lambda=linspace(1,5000,2000)*1000;

[LL,RR]=meshgrid(lambda,r);
KK4=(2*pi./LL).^4;

% Forsyth Eqs. (3) and (6)
xai=1+D.*KK4/drho/g;
phi=1+D.*KK4/rc/g;
upco=exp(-(2*pi./LL)*zm);

% Olhede and Simons Eq. 59
fax=sqrt(f2)*rc./drho;

% Combined admittance for independent surface
% and subsurface loading
% Olhede and Simons Eq. (60), same thing
% Could write this using ADMITTANCE
QTBRP=-2*pi*G*rc*upco.*(xai+fax.^2.*phi-RR.*fax.*(phi.*xai+1))...
     ./(xai.^2+fax.^2-2*RR*fax.*xai);
% Note this has a root at:
lsing=2*pi/[(max(real(roots([1 -2*fax fax^2])))-1)*g*drho./D]^(1/4);

% Now for the negative correlation coefficients 
% nneg=find(~sum(RR==0,2));
RR=-RR;
QTBRN=-2*pi*G*rc*upco.*(xai+fax.^2.*phi-RR.*fax.*(phi.*xai+1))...
     ./(xai.^2+fax.^2-2*RR*fax.*xai);
% and also -1
QTBRN=[QTBRN ; -2*pi*G*rc*upco(1,:).*(xai(1,:)+...
      fax.^2.*phi(1,:)-(-1).*fax.*(phi(1,:).*xai(1,:)+1))...
     ./(xai(1,:).^2+fax.^2-2*(-1)*fax.*xai(1,:))];

% Convert to mgal/m
QTBRP=QTBRP/1e-5;
QTBRN=QTBRN/1e-5;

% Cross check with ADMITTANCE which itself has been cross-checked with FORSYTH
% for uncorrelated loading only
tzer=find(sum(RR==0,2));
difer(QTBRP(tzer,:)-admittance(Te,f2,0,zm,rc,drho,[],2*pi./lambda,0),...
      [],[],NaN);
difer(QTBRN(tzer,:)-admittance(Te,f2,0,zm,rc,drho,[],2*pi./lambda,0),...
      [],[],NaN);
% Now do the full cross-check with ADMITTANCE
for  in=1:length(r)
  difer(QTBRP(in,:)-admittance(Te,f2,r(in),zm,rc,drho,[],2*pi./lambda,0),...
      [],[],NaN);
  difer(QTBRN(in,:)-admittance(Te,f2,-r(in),zm,rc,drho,[],2*pi./lambda,0),...
      [],[],NaN);
end
% Test -1
difer(QTBRN(in+1,:)-admittance(Te,f2,-1,zm,rc,drho,[],2*pi./lambda,0),...
      [],[],NaN);

% Do the plotting proper
figure(gcf)
for index=1:length(r)
  ptp(index)=semilogx(lambda/1000,QTBRP(index,:),'k-');  
  hold on
end
for index=1:length(r)+1
  ptn(index)=semilogx(lambda/1000,QTBRN(index,:),'k--');
end
hold off

% Cosmetics
set(gca,'xdir','rev','ydir','rev')
yls=[-0.2 0.01]+0.05;
ylim(yls)
xlim([10 5000])
hold on
ppp=semilogx([lsing lsing]/1000,ylim,'--','linew',2);
set(ppp,'Color',grey)
hold off
bottom(ppp,gca)

grid on
longticks

%xl=xlabel('Wavelength (km)');
%yl=ylabel('Bouguer per Topography: Admittance (mgal/m)');
xl=xlabel('wavelength (km)');
yl=ylabel('admittance Q_\circ (mgal/m)');

set(ptp,'linew',2)
set(ptn,'linew',2)

% Legend boxes
pos=[0.5100    0.5100
     0.4950    0.4600
     0.4800    0.4100
     0.4650    0.3600
     0.4650    0.3100
     0.4500    0.2600
     0.4350    0.2100
     0.4200    0.1600
     0.0400    0.0465];

hold on
b=gca;
d=axes;
tbb=0.015*infl;
tbx=0.03*onfl;
fb=fillbox(...
    [pos(1:end-1,1)-tbx pos(1:end-1,1)+tbx pos(1:end-1,2)+tbb pos(1:end-1,2)-tbb],'w');
hold on
fb=[fb ; fillbox(...
    [pos(end,1)-tbx/onfl pos(end,1)+tbx/onfl pos(end,2)+tbb pos(end,2)-tbb],'w')];

axes(d)
axis([0 1 0 1])
set(d,'color','none')
nolabels(d)
noticks(d)

set(fb(1:length(r)),'facecolor',[0.75 0.75 0.75])

% Now plot the text boxes from the top, most negative down
ers=unique([r -r -1]);
for index=1:length(ers)
  t(index)=text(pos(index,1),pos(index,2),num2str(ers(index)));
end
t(index+1)=text(pos(index+1,1),pos(index+1,2),'r');

set(t,'horizontala','center')

set(d,'position',get(b,'position'))

set(b,'XTickLabel',[10 100 1000])
set([gca xl yl],'fonts',12)

defval('land',1)
if land==1
  fig2print(gcf,'landscape')
  figdisp
end

varns={d,xl,yl,ptp,ptn,t};
varargout=varns(1:nargout);

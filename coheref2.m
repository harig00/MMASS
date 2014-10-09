function varargout=coheref2(infl,onfl,land)
% COHEREF2(infl,onfl,land)
%
% Expected coherence from Forsyth (1985)
% One elastic thickness, several loading ratios.
%
% Forsyth (1985), Figure 8.
%
% OUTPUT:
%
% Lots of handle variables
%
% See also: ADMITR, ADMITB, ADMITF2, COHERETE
%
% Last modified by fjsimons-at-alum.mit.edu, 03/07/2012

defval('infl',1)
defval('onfl',1)

% Elastic thickness [m]
Te=40*1000;

% Loading ratios
f2=[0.1 0.5 1 2 5];

% Initial-loading correlation coefficient
r=0;

% Wavelengths [m]
lambda=linspace(10,2000,1000)*1000;

% Do the calculations
C=forsyth(Te,lambda,f2,r);

% Make the plot, wavelength in km
p=semilogx(lambda/1000,C,'k');
p=p(:)';

set(gca,'xdir','rev')

yls=[-0.025 1.025];
ylim(yls)
xlim([10 2000])

%xl(1)=xlabel('Wavelength (km)');
%yl(1)=ylabel('Bouguer Anomaly/Topography Coherence \gamma^2');
xl(1)=xlabel('wavelength (km)');
yl(1)=ylabel('coherence \gamma_f^2');

grid on

set(p,'linew',2)

post=[0.3200    0.40  % 0.1
      0.3300    0.45  % 0.5
      0.3400    0.50  % 1
      0.3500    0.55  % 2
      0.3600    0.60  % 5
      0.0400    0.04];

post(1:5,1)=post(1:5,1)-0.05;
  
hold on
b=gca;
d=axes;
tbb=0.015*infl;
tbx=0.03*onfl;
fb=fillbox(...
    [post(1:end-1,1)-tbx post(1:end-1,1)+tbx post(1:end-1,2)+tbb post(1:end-1,2)-tbb],'w');
hold on
fb=[fb ; fillbox(...
    [post(end,1)-tbx/onfl post(end,1)+tbx/onfl post(end,2)+tbb*1.25 post(end,2)-tbb*1.25],'w')];
axes(d)
axis([0 1 0 1])
set(d,'color','none')
noticks(d)
nolabels(d)

for index=1:length(f2)
  t(index)=text(post(index,1),post(index,2),num2str(f2(index)));
end
t(index+1)=text(post(index+1,1),post(index+1,2),'f^2');

set(t,'horizontala','center')

set(d,'position',get(b,'position'))

longticks(b)

set(b,'XTickLabel',[10 100 1000])

set([gca xl yl],'fonts',12)

defval('land',1)
if land==1
  fig2print(gcf,'landscape')
  figdisp
end

varns={d,xl,yl,p,t};
varargout=varns(1:nargout);



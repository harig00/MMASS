function varargout=coherete(infl,onfl,land)
% COHERETE(infl,onfl,land)
%
% Expected coherence from Forsyth (1985)
% One loading ratio, several elastic thicknesses.
%
% Forsyth (1985), Figure 7.
%
% OUTPUT:
%
% Lots of handle variables
%
% See also: ADMITR, ADMITB, ADMITF2, COHEREF2
%
% Last modified by fjsimons-at-alum.mit.edu, 03/07/2012

defval('infl',1)
defval('onfl',1)

% Elastic thicknesses [m]
Te=[1 10:10:100]*1000;

% Loading ratio
f2=[1];

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

% xl(1)=xlabel('Wavelength (km)');
% yl(1)=ylabel('Bouguer Anomaly/Topography Coherence \gamma^2');
xl(1)=xlabel('wavelength (km)');
yl(1)=ylabel('coherence \gamma_f^2');

grid on

set(p,'linew',2)

post=[0.1814    0.3000
      0.1904    0.35
      0.2016    0.4000
      0.2161    0.45
      0.2307    0.5
      0.2520    0.55
      0.2766    0.6
      0.3102    0.65
      0.3617    0.70
      0.4535    0.7
      0.7906    0.7000
      0.040     0.04];
  
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

for index=1:length(Te)
  t(index)=text(post(index,1),post(index,2),...
		num2str(Te(length(Te)-index+1)/1000));
end
t(index+1)=text(post(index+1,1),post(index+1,2),'Te');

set(t,'horizontala','center')

longticks(b)

set(b,'XTickLabel',[10 100 1000])

set([gca xl yl],'fonts',12)

set(d,'position',get(b,'position'))

defval('land',1)
if land==1
  fig2print(gcf,'landscape')
  figdisp
end

varns={d,xl,yl,p,t};
varargout=varns(1:nargout);


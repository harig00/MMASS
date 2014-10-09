function varargout=coherer(infl,onfl,land)
% COHERER(infl,onfl,land)
%
% Expected coherence from Forsyth (1985)
% One elastic thickness, several loading correlations.
%
% See Crosby (2007)
%
% OUTPUT:
%
% Lots of handle variables
%
% See also: ADMITR, ADMITB, ADMITF2, COHERETE, COHEREF2
%
% Last modified by fjsimons-at-alum.mit.edu, 03/07/2012

defval('infl',1)
defval('onfl',1)

% Elastic thickness [m]
Te=40*1000;

% Loading ratio
f2=1;

% Correlation coefficients (r and -r will be considered)
r=[0 0.25 0.5 0.75];

% Wavelengths [m]
lambda=linspace(10,2000,1000)*1000;

% Do the calculations
CP=forsyth(Te,lambda,f2,r);
CM=forsyth(Te,lambda,f2,-r);

% Check that +/-1 is not worth showing
difer(1-forsyth(Te,lambda,f2,-1),[],[],NaN);
difer(1-forsyth(Te,lambda,f2,+1),[],[],NaN);

% Make the plot, wavelength in km
pp=semilogx(lambda/1000,CP,'k');
pp=pp(:)';
hold on
pm=semilogx(lambda/1000,CM,'k--');
pm=pm(:)';

set(gca,'xdir','rev')

yls=[-0.025 1.025];
ylim(yls)
xlim([10 2000])

%xl(1)=xlabel('Wavelength (km)');
%yl(1)=ylabel('Bouguer Anomaly/Topography Coherence \gamma^2');
xl(1)=xlabel('wavelength (km)');
yl(1)=ylabel('coherence \gamma_\circ^2');

grid on

set([pp pm],'linew',2)
post=[0.3531    0.6932 % -0.75
      0.3040    0.6181 % -0.5
      0.2903    0.5376 % -0.25
      0.2848    0.4445 %  0
      0.2739    0.3658 %  0.25
      0.2671    0.2657 %  0.5
      0.2630    0.1351 %  0.75
      0.0400    0.04];

hold on
b=gca;
d=axes;
tbb=0.015*infl;
tbx=0.03*onfl;
fb=fillbox(...
    [post(1:end-1,1)-tbx post(1:end-1,1)+tbx post(1:end-1,2)+tbb post(1:end-1,2)-tbb],'w');
hold on
fb=[fb ; fillbox(...
    [post(end,1)-tbx/onfl post(end,1)+tbx/onfl post(end,2)+tbb post(end,2)-tbb],'w')];
axis([0 1 0 1])
axes(d)
set(d,'color','none')
noticks(d)

set(fb(1:length(r)-1),'facecolor',[0.75 0.75 0.75])

ers=unique([r -r]);
for index=1:length(ers)
  t(index)=text(post(index,1),post(index,2),num2str(ers(index)));
end
t(index+1)=text(post(index+1,1),post(index+1,2),'r');

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

varns={d,xl,yl,pp,pm,t};
varargout=varns(1:nargout);



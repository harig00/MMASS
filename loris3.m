function loris3(opt,agu,tperc)
% LORIS3(opt,agu,tperc)
%
% Plots a wavelet based multiscale rendition of topography for a single
% face of the cubed sphere using different orthonormal wavelets. Also
% plots histograms of the wavelet coefficients, and the like.
%
% INPUT:
%
% opt    0 plots wavelet coefficients and histograms
%        1 plots map reconstructions and histograms
%        2 plots wavelet coefficients, maps, and histograms [default]
% agu    1 presentation-type plot of a single transform, only operational
%        when opt is the default 2
%
% Last modified by fjsimons-at-alum.mit.edu, 4/10/2012

defval('opt',2)
defval('agu',0)
% Truncation level as percentiles, passed on to LORISTOPO2D
% For the moment only to D4
defval('tperc',85);
% For SPIE this was 
defval('iface',1);
% For GJI this was 3
defval('iface',5);

clf
if opt<2
  [ah,ha,H]=krijetem(subnum(2,3));
else
  [ah,ha,H]=krijetem(subnum(3,3));
end

if agu==1
  aho=ah;
  ah=ha;
  ha=aho;
end

% The extra percentiles for the histogram plot
experc=[0.075 100];
% The number of bins in the histograms
numbin=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(1))
if opt<2
  [hd2,fd2e,gd2,jd2,tt2,rd2,t2,vw2]=loristopo2d('D2',opt,[],[],[],[],iface);
else
  [hd2,fd2e,gd2,jd2,tt2,rd2,t2,vw2]=loristopo2d('D2',1,[],[],[],[],iface);
  delete(tt2)
end

% These are the percentiles you will be getting
percs2=unique([experc t2(2,:)]);

axes(ha(2))
if opt==2
  [hd2,fd2e,gd2b,jd2,tt2,rd2,t2,vw2]=loristopo2d('D2',2,[],[],[],[],iface);
  axes(ha(3))
end
[hi2pl,hi2pr,hi2ml,hi2mr,yr2,pv2]=...
    newhist('D2',vw2,numbin,rd2(3),percs2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))
if opt<2
  [hd4,fd4,gd4,jd4,tt4,rd4,t4,vw4]=loristopo2d('D4',opt,[],[],[],[],iface);
else
  [hd4,fd4,gd4,jd4,tt4,rd4,t4,vw4]=...
      loristopo2d('D4',1,[],[],[],[],iface,[],tperc);
  delete(tt4)
end

% These are the percentiles you will be getting
percs4=unique([experc t4(2,:)]);

axes(ha(4))
if opt==2
  axes(ha(5))
  [hd4,fd4,gd4b,jd4b,tt4,rd4,t4,vw4]=...
      loristopo2d('D4',2,[],[],[],[],iface,[],tperc);
  axes(ha(6))
end
[hi4pl,hi4pr,hi4ml,hi4mr,yr4,pv4]=...
    newhist('D4',vw4,numbin,rd4(3),percs4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(3))
if opt<2
  [hd6,fd6,gd6,jd6,tt6,rd6,t6,vw6]=loristopo2d('D6',opt,[],[],[],[],iface);
else
  [hd6,fd6,gd6,jd6,tt6,rd6,t6,vw6]=loristopo2d('D6',1,[],[],[],[],iface);
  delete(tt6)
end

% These are the percentiles you will be getting
percs6=unique([experc t6(2,:)]);

axes(ha(6))
if opt==2
  axes(ha(8))
  [hd6,fd6,gd6b,jd6,tt6,rd6,t6,vw6]=loristopo2d('D6',2,[],[],[],[],iface);
  axes(ha(9))
end
[hi6pl,hi6pr,hi6ml,hi6mr,yr6,pv6]=...
    newhist('D6',vw6,numbin,rd6(3),percs6);

% Cosmetics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
longticks(ah(length(ah)-2:length(ah)))
groy=1.05;
set(ah(length(ah)-2:length(ah)),'ylim',...
		  [-max([yr2 yr4 yr6]) max([yr2 yr4 yr6])]*groy)
grox=50;
xlim2=minmax(pv2)+[0 range(pv2(:))/grox];
xlim4=minmax(pv4)+[0 range(pv4(:))/grox];
xlim6=minmax(pv6)+[0 range(pv6(:))/grox];
set(ah(length(ah)-2),'xlim',xlim2)
set(ah(length(ah)-1),'xlim',xlim4)
set(ah(length(ah)-0),'xlim',xlim6)
% Set the tick marks to the 5, 50 and 95 percentiles... do check as I am
% putting in their location in an absolute sense for my convenience. As
% the positive and negative parts have a very similar distribution, we
% take their average. This too, should be checked for truthfulness.
wper=2:4;
display(sprintf('Histogram tick marks at average %i %i %i percentiles',...
		percs2(wper)))
display(sprintf('Histogram tick marks at average %i %i %i percentiles',...
		percs4(wper)))
display(sprintf('Histogram tick marks at average %i %i %i percentiles',...
		percs6(wper)))
tix2=unique([xlim2 mean(pv2(:,wper))]);
tix4=unique([xlim4 mean(pv4(:,wper))]);
tix6=unique([xlim6 mean(pv6(:,wper))]);
% And include the end points of the current axis also
set(ah(length(ah)-2),'xtick',tix2,'xtickl',round(10*tix2)/10,'xgrid','on')
set(ah(length(ah)-1),'xtick',tix4,'xtickl',round(10*tix4)/10,'xgrid','on')
set(ah(length(ah)-0),'xtick',tix6,'xtickl',round(10*tix6)/10,'xgrid','on')

if opt<2
  sr=0.785;
  shrink(ah([1 2 3]),sr,sr)
  movev(ah([1 2 3]),-0.05)
  movev([gd2 gd4 gd6],-0.075)
  shrink([gd2 gd4 gd6],sr,1.5)
  shrink(ah([4 5 6]),1,1.2)
  movev(ah([4 5 6]),.05)
  fig2print(gcf,'portrait')
else
  if agu~=1
    sr=1;
    shrink([gd2 gd4 gd6],sr,1.5)
    shrink([gd2b gd4b gd6b],sr,1.5)
    shrink(ah(length(ah)-2:length(ah)),1,1.475)
    shrink(ah(length(ah)-2:length(ah)),1.1,1)
    movev(ah([1 2 3]),-0.025)
    movev([gd2 gd4 gd6],-0.025)
    movev(ah(length(ah)-2:length(ah)),0.075)
    fig2print(gcf,'tall')
  elseif agu==1
    fig2print(gcf,'landscape')
  end
end

if opt==2
  % Move them close together
  movev(tt2,7.5)
  movev(tt4,7.5)
  movev(tt6,7.5)
  movev(ah(7:9),0.01)
end

figdisp([],opt,[],0)

if agu==1
  delete(ha(1:3))
  delete(ha(7:9))
  delete([gd2 jd2 gd6 jd6 gd2b gd6b])
  sr=0.5;
  shrink(ah(2),sr,sr)
  shrink([gd4 ],sr,1.25)
  movev(gd4,-.075)
  movev(tt4,-10)
  shrink(ha(5),sr,sr)
  shrink([gd4b],sr,1.25)
  movev(gd4b,-.075)
  shrink(ha(6),1/0.91,0.85)
  movev(ha(6),.0325)
  figdisp([],sprintf('%i_%i_%i',opt,agu,round(tperc)),[],0)
  moveh([ah(2) gd4],.05)
  moveh([ha(5) gd4b],.025)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hipl,hipr,himl,himr,yrg,percval]=...
       newhist(wav,data,numbin,trunx,percs)
% This gives us logarithmic histograms of positive and negative values
% separately. The bins terminate at the truncation so we can color them

defval('percval',[5 50 95])

N=prod(size(data));

% Divide the data into signed ranges
dp=log10(data(data>0));
dm=log10(abs(data(data<0)));
disp(sprintf('%i are exactly zero and are lost',sum(data(:)==0)))

% And here we've been truncating hard
dt=log10(trunx);

% This will do this without regards of the truncation point
% [ap,bp]=hist(dp,numbin);
% hip=bar(bp,ap,1); hold on
% [am,bm]=hist(dm,numbin);
% him=bar(bm,-am,1); hold off

% This will do this while respecting the truncation point
[binpl,binpr]=newbins(dp,dt,numbin);
% Their sum equals the number of positive values
apl=histc(dp,binpl);
apr=histc(dp,binpr);
hipl=bar(binpl,apl/N*100,'histc'); hold on
hipr=bar(binpr,apr/N*100,'histc'); 

[binml,binmr]=newbins(dm,dt,numbin);
% Their sum equals the number of negative values
aml=histc(dm,binml);
amr=histc(dm,binmr);
himl=bar(binml,-aml/N*100,'histc'); 
himr=bar(binmr,-amr/N*100,'histc'); hold off

% Check that the total sum of everything adds up to a hundred percent 
difer(sum([aml; amr; apl; apr]/N*100)-100,[],[],NaN)

% Cosmetics
delete(findobj('marker','*'))
set(himl,'FaceC','w')
set(himr,'FaceC','r')
set(hipl,'FaceC','w')
set(hipr,'FaceC','b')
xlim(minmax([dp ; dm]))

% Collect the y-range
yrg=minmax([apl ; apr; aml ; amr]/N*100);

% Collect percentiles for output only
% Later one we will mark the mean
disp(sprintf(sprintf('Calculating %s percentiles',...
		     repmat('%5.3f ',1,length(percs))),percs))
percval=[prctile(dp,percs) ; prctile(dm,percs); prctile([dm ; dp],percs)];

strunk=sprintf('log_{10}(%s coefficients)',wav);
xl=xlabel(strunk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [binl,binr]=newbins(d,dt,numbin)
numl=ceil([dt-min(d)]/range(d)*numbin);
numr=floor([max(d)-dt]/range(d)*numbin);
binl=linspace(min(d),dt,numl+1);
binr=linspace(dt,max(d),numr+1);

function dsbpboxcap
% DSBPBOXCAP
%
% Makes a plot of the boxcar power spectrum
% Simons & Dahlen (2007), Figure 2
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2006

% For some reason, the second time round, the plot looks right, but not
% the first... rerun!

TH1=[10 20 30];
[Bl1,dels]=bpboxcap(TH1,100,[],0,1);

TH2=[30 20 10];
[Bl2,dels]=bpboxcap(TH2,100,[],0,2);
Bl2(2:2:end,:)=NaN;

for index=1:length(TH1)
  % Note that this is COLUMN per COLUMN !!
  Bl1(:,index)=decibel(Bl1(:,index));
  Bl2(:,index)=decibel(Bl2(:,index));
  % Make sure the max is the B0 term
  difer(max(Bl1(:,index))-Bl1(1,index))
  difer(max(Bl2(:,index))-Bl2(1,index))
end

% Make the figure
clf
[ah,ha]=krijetem(subnum(3,2));
minb=min([Bl1(:) ; Bl2(:)]);
maxb=-minb;
tix=[0:20:maxb];
widr=[0 5];

% Single cap
for index=1:3
  axes(ha(index))
  b1(index)=bar(dels,Bl1(:,index)-minb,1);
  set(gca,'ytick',sort(maxb-tix),'ytickl',fliplr(-tix))
  set(gca,'ylim',minmax(Bl1(:,index)-minb)+widr)
  xl(index)=xlabel('degree p');
  yl(index)=ylabel('power B_p/B_0 (dB)');
end

% Double cap
for index=1:3
  axes(ha(index+3))
  b2(index)=bar(dels,Bl2(:,index)-minb,1);
  set(gca,'ytick',sort(maxb-tix),'ytickl',fliplr(-tix))
  set(gca,'ylim',minmax(Bl2(:,index)-minb)+widr)
  xl(index+3)=xlabel('degree p');
  yl(index+3)=ylabel('power B_p/B_0 (dB)');
end

% Cosmetics
set(ah,'xlim',[-1 101],'ylim',[30 maxb]+widr)
set([b1 b2],'FaceC',grey)
longticks(ah)
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ah(5:6),1/2,'across')
% For some reason this doesn't work
serre(ha(1:3),1/3,'down')
serre(ha(4:6),1/3,'down')
% Leave to last
nolabels(ah(1:4),1)
nolabels(ha(4:6),2)
delete(xl([1 2 4 5]))
delete(yl([4:6]))

% Make miraculous adjustment
set(ah(6),'position',...
	  [getpos(ah(6),1) getpos(ah(5),2) ...
	   getpos(ah(6),3) getpos(ah(5),4)])

% Logo scale and location
logsk=20;
cloc='ll';

% Put in equivalent wavelength labels
for index=1:3
  % Get equivalent wavelengths for the single cap
  for N=1:5
    skel(N)=max(roots([1 1 -(N*180./TH1(index))^2]));
  end
  % Put in cap size labels
  legsi{index}=sprintf('  %s = %i%s','\Theta',TH1(index),str2mat(176));
  [bh(index),th(index)]=boxtex('ur',ha(index),legsi{index},12);
  % Put on top labels
  [ax(index),axl(index)]=...
      xtraxis(ha(index),skel,1:5,'number of wavelengths in \itR');
  % Put the cap logo on
  lah(index)=caplogo(ha(index),1,cloc,1/logsk,1/logsk,20);
end
clear skel
EN=1:2:5;
for index=1:3
  % Get equivalent wavelengths for the double cap
  for N=1:length(EN)
    skel(N)=max(roots([1 1 -(EN(N)*180./(90-TH2(index)))^2]));
  end
  % Put in cap size labels
  legsi{index+3}=sprintf('  %s = %i%s','\Theta',90-TH2(index),str2mat(176));
  [bh(index+3),th(index+3)]=boxtex('ur',ha(index+3),legsi{index+3},12);
  % Put on top labels
  [ax(index+3),axl(index+3)]=...
      xtraxis(ha(index+3),skel,EN,'number of wavelengths in \itR');
  % Put the cap logo on
  [lah(index+3),xal(index+3)]=caplogo(ha(index+3),3,cloc,1/logsk,1/logsk,...
				      20);
end
set(ax,'xgrid','on')
set(ah,'ygrid','on')
delete(axl([2 3 5 6]))
longticks(ax)

% Make the extra axis come last
for index=1:length(ah)
  axes(ah(index))
  set(ah(index),'color','none')
end
for index=1:length(lah)
  axes(lah(index))
end

set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')
figdisp

disp('Rerun this for the correct label dimensions!')

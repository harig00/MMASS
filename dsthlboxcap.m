function dsthlboxcap(grayscale)
% DSTHLBOXCAP(grayscale)
%
% Makes surface plots of the power spectrum of the boxcar caps.
% Dahlen & Simons (2008), Figure 3
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2006

% For some reason, the second time round, the plot looks right, but not
% the first... rerun!

defval('grayscale',64)

defval('THint',0.1)

% Later on, for double cap, do fliplr instead of 90-TH
TH=THint:THint:90-THint;

% Get the boxcar power spectrum
[Bl1,dels]=bpboxcap(TH,100,[],0,1);
[Bl2,dels]=bpboxcap(TH,100,[],0,2);

% COLUMN per COLUMN in decibel is what we want
warning off
for index=1:length(TH)
  % Note that this is COLUMN per COLUMN !!
  Bl1(:,index)=decibel(Bl1(:,index));
  Bl2(:,index)=decibel(Bl2(:,index));
end
% Make sure the max was indeed the B0 term
difer(max(Bl1)-Bl1(1,:))
difer(max(Bl2)-Bl2(1,:))
% Make the odd degrees invisible rather than -Inf
Bl2(2:2:end,:)=NaN;
warning on

% Take out the tiny values at 90 or else...
cax1=[-60 0];

% Commence plotting
clf
ah=krijetem(subnum(1,2));

axes(ah(1))
imagefnan([TH(1) dels(1)],[TH(end) dels(end)],Bl1,...
	  flipud(gray(grayscale)),cax1,[],[],0)
axis ij
xlim([0 90])
axis square
yl(1)=ylabel('degree p');
xl(1)=xlabel(sprintf('single cap radius %s','\Theta'));
tl(1)=title('boxcar power B_p/B_0 (dB)');

% Maybe hold on and plot the lines on?
warning off
hold on
for index=[1 2 3 4 5]
  thel=index*180./sqrt(dels.*(dels+1));
  thel(1)=NaN;
  pb(index)=plot(thel,dels,'color','w','linew',0.5);
end
warning on

% Plot equivalent wavelengths
coords2=[ 9   17
	 14   22
	 19   27
	 23   32
	 26   37];
hold on
f1=10;
f2=8;
for index=1:size(coords2,1)
  [b2(index),t2(index)]=boxtex([coords2(index,1) coords2(index,2)],gca,...
			       sprintf('%i %s',index,'\times'),f1,0.7,0.6);
end
set(t2,'FontS',f2)
set(t2,'Color','k')
set(b2,'EdgeC','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))
imagefnan([TH(1) dels(1)],[TH(end) dels(end)],fliplr(Bl2),...
	  flipud(gray(grayscale)),cax1,[],[],0)
axis ij
xlim([0 90])
axis square
yl(2)=ylabel('degree p');
xl(2)=xlabel(sprintf('double cap radius %s','\Theta'));
tl(2)=title('boxcar power B_p/B_0 (dB)');
warning off
hold on
for index=[1 2 3 4 5]
  thel=index*180./sqrt(dels.*(dels+1));
  thel(1)=NaN;
  pb2(index)=plot(thel,dels,'color','k','linew',0.5);
end
warning on
% Plot equivalent wavelengths
coords2=[ 9   17
	 14   22
	 19   27
	 23   32
	 26   37];
hold on
f1=10;
f2=8;
for index=1:size(coords2,1)
  [b2(index),t2(index)]=boxtex([coords2(index,1) coords2(index,2)],gca,...
			       sprintf('%i %s',index,'\times'),f1,0.7,0.6);
end
set(t2,'FontS',f2)
set(t2,'Color','k')
set(b2,'EdgeC','k')

% Cosmetics
fig2print(gcf,'portrait')

% Color bar
colormap(flipud(gray(grayscale)))
c=colorbar('hor');
shrink(c,1,2)
movev(c,-0.1)
moveh(c,-0.2)
longticks(c,2)
clims=get(c,'xlim');
set(c,'xtick',linspace(clims(1),clims(end),max(abs(cax1))/10+1),...
      'xtickl',cax1(1):10:cax1(end))
axes(c)
xcl=xlabel('boxcar power B_p/B_0 (dB)');
movev(xcl,-1.25)

% Cosmetics
longticks(ah)
delete(yl(2))
nolabels(ah(2),2)
set(ah,'xtick',[0:15:90])
serre(ah(1:2),1/2,'across')
deggies(ah,1)

% Put the cartoon on...
% Plot the logo on it
l(1)=caplogo(ah(1),1,[],[],[],20);
shrink(l(1),1/2,1/2)
movev(l(1),0.21)
moveh(l(1),-0.015)

l(2)=caplogo(ah(2),3,[],[],[],20);
shrink(l(2),1/2,1/2)
movev(l(2),0.21)
moveh(l(2),-0.015)

delete(tl)

set(gcf,'color','w','inverthardcopy','off')
figdisp


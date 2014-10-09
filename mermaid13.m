function mermaid13
% MERMAID13
%
% Last (maybe not) MERMAID paper figure - zooms on the detect onsets
% A mangled version of MERMAID03. Also with wavelet denoising...
%
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2008

dirs={'11042003','11042003','11042003',...
      '09102004','08092007'};
files={'mpilot1_21118','mpilot1_98947','mpilot1_138106',...
       'mpilot2_79891','mpilot3_4759'};
filts={'lowpass','bandpass','bandpass','bandpass','bandpass'};
fops={{1,2,2},{1,8,2,2},{1,8,2,2},{1,8,2,2},{1,8,2,2}};

% Some arbitrary y-axis scalings
scax=[1e3 1e4 1e4 1e3 5e4];
scay=[1 1 1.2 1 1.5];
scaz=[1e3 1e4 1e4 1e3 5e4];

clf
[ah,ha]=krijetem(subnum(5,4));

% New in MERMAID13
xadj=[-5 -13.4 -15 +20 +4];
xtox=[-15:5:15]; % Symmetric!

for inx=1:5
  disp(' ')
  dir=fullfile('/home/fjsimons/MERMAID',dirs{inx},'EVENTS');
  
  % Time-domain plots of the identified events %%%%%%%%%%%%%%%%%%%%
  axes(ha(inx))
  [x,h,Fs,p(inx),xl(inx),yl(inx),tl(inx),pt01(inx,:),xf]=...
      timdomplot(fullfile(dir,sprintf('%s.sac',files{inx})),2,...
		 filts{inx},fops{inx});
  % Adjust axes to be pretty - NOW THIS IS DIFFERENT FROM MERMAID03
  % And with this zoom you notice pretty serious misalignment - so, realign
  xtix=xtox+xadj(inx);

  txl=xtix(1)+range(xtix)/60;
  txr=xtix(end)-range(xtix)/60;

  set(p(inx),'ydata',get(p(inx),'ydata')/scax(inx))
  set(pt01(inx,1),'ydata',get(pt01(inx,1),'ydata')/scax(inx))
  set(pt01(inx,2),'ydata',get(pt01(inx,2),'ydata')/scax(inx))
  set(ha(inx),'ylim',get(pt01(inx,2),'ydata')/scay(inx))
  % Now the zero is the eyeballed, not the predicted, P wave arrival
  set(ha(inx),'xtick',h.T0+xtix,'xtickl',xtix-xadj(inx),'xlim',h.T0+[xtix(1) xtix(end)])
  bottom(pt01(inx,1),ha(inx))
  bottom(pt01(inx,2),ha(inx))
  % Now add magnitude, event depth and distance information on these plots
  ylo=ylim;
  leg1=sprintf('%s = %3.1f%s','\Delta',h.GCARC,str2mat(176));
  leg2=sprintf('mag = %3.1f ',h.MAG);
  % And the filter "compression" legend
  leg3=sprintf('%i bps',16*round(1/h.DELTA));
  
  tx1(inx)=text(h.T0+txr,ylo(1)+range(ylo)/10,leg1,...
		'horizontala','right');
  tx2(inx)=text(h.T0+txr,ylo(2)-range(ylo)/10,leg2,...
		'horizontala','right');
  tx3(inx)=text(h.T0+txl,ylo(2)-range(ylo)/10,leg3,...
		'horizontala','left');

  % Time-domain plots of WAVELET denoised signals %%%%%%%%%%%%%%%%%%%%
  % Do the wavelet analysis - on the raw data x or on filtered xf?
  % Copy the following from SCALOGRAMPLOT
  defval('tipo','CDF')
  defval('nvm',[2 4])
  defval('nd',5)
  defval('pph',3)
  defval('intel',0)

  % What are we wavelet transforming? The raw signal or the filtered signal?
  xw=xf; % DEFINITELY the filtered signal!!

  % Actually, let's see to it here that we need to only take a certain
  % portion of it!!
  allofx=linspace(h.B,h.E,h.NPTS);
  portionx=h.T0+[xtix(1) xtix(end)];
  
  % Segment only what will be transformed
  xw=xw(allofx>=portionx(1) & allofx<=portionx(2));
  partofx=allofx(allofx>=portionx(1) & allofx<=portionx(2));

  % Now also figure out the signal-to-noise ratio of the first "half" 
  % compared to the second - possibly the halves are not symmetric
  % Assume symmetry, or rather check for it
  if xtox(1)==-xtox(end)
    % plot(1:length(xr)/2,xr(1:length(xr)/2)); hold on
    % SNR of the Original portion
    dblev0=decibel(mean(xw(1+length(xw)/2:end).^2),...
		   mean(xw(1:length(xw)/2).^2));
  end

  legD=sprintf('SNR = %i dB',round(dblev0));
  txD(inx)=text(h.T0+txl+0,ylo(1)+range(ylo)/10,legD,...
		'horizontala','left');

  axes(ha(inx+5))

  ochk=1*mod(length(xw),2);
  % Don't for
  [a,d,an,dn]=wt(detrend(xw(1:end-ochk)),tipo,nvm,nd,pph,intel);
  oldsum=an(end)+sum(dn); % Remember polyphase lengthens slightly
  % The division factor is pretty crucial here, and use across all scales
  % [d,dn]=threshold(d,dn,'soft',3);
  [d,dn]=threshold2(d,dn,'soft',1);
  newsum=an(end)+sum(dn);
  [xc,xr,ts]=iwt(a,d,an,dn,tipo,nvm,pph);
  % cf=compres(x,xc,an,dn);
  cf=newsum/length(xw)*100;

  pxr(inx)=plot(partofx,xr/scaz(inx),'k');

  % Same bag of tricks as in the column on the left
  set(ha(inx+5),'xtick',h.T0+xtix,'xtickl',xtix-xadj(inx),'xlim',portionx)
  % Same axis limits left and right
  set(ha(inx+5),'ylim',get(ha(inx),'ylim'))

  % I want to quote the COMPRESSED rate but only for the displayed segment
  leg4=sprintf('%i bps (%i %s)',round(16*round(1/h.DELTA)*cf/100),round(cf),'%');
  
  tx4(inx)=text(h.T0+txl,ylo(2)-range(ylo)/10,leg4,...
		'horizontala','left');

  xl(inx+5)=xlabel('time (s)');

  % Now also figure out the signal-to-noise ratio of the first "half" 
  % compared to the second - possibly the halves are not symmetric
  % Assume symmetry, or rather check for it
  if xtox(1)==-xtox(end)
    % plot(1+length(xr)/2:length(xr),xr(1+length(xr)/2:end),'r')
    % SNR of the Wavelet-denoised portion of it
    dblevW=decibel(mean(xr(1+length(xr)/2:end).^2),...
		   mean(xr(1:length(xr)/2).^2));
  end

  legD=sprintf('SNR = %i dB',round(dblevW));
  txD(inx+5)=text(h.T0+txl+0,ylo(1)+range(ylo)/10,legD,...
		'horizontala','left');
end

% Move everything up to make room for the colorbars
movev(ha,.1)
set(p,'Color','k')
delete(pt01)
for inx=1:2
  serre(ha([1:5]+(inx-1)*5),1/2,'down')
end
serre(ah([1 2 ; 5 6 ; 9 10; 13 14; 17 18]))

delete(ha(11:end))
ah=ah([1 2 5 6 9 10 13 14 17 18]);
ha=ha(1:10);

% Cosmetic adjustments
% Mark the P-wave arrival time
for inx=1:2
  % Here need to make sure that there are as many blanks as needed
  axx(inx)=xtraxis(ha(1+(inx-1)*5),get(ha(1+(inx-1)*5),'xtick'),...
	      {'' '' '' 'P' '' '' ''});
end
longticks([ah axx])

set([ha xl yl axx tx1 tx2 tx3 tx4 txD],'FontS',8)

set([p pxr],'LineW',0.5)
delete(xl([1:4 6:9]))
set(xl(5),'string','time (s)')
set(yl,'string','sound pressure (scaled hydrophone counts)')
for inx=1:2
  nolabels(ha([1:4]+(inx-1)*5),1) % This must be first
end
nolabels(ha(6:10),2)
delete(tl)
[bh,th]=label(ha,'ll',10,[],[],[],[],0.75);

% Now get the extent of the added labels and shift legD by that amount
for inx=1:10
  shifit=diff(unique(get(bh(inx),'XData')));
  set(txD(inx),'Position',getpos(txD(inx))+[1.2*shifit 0 0])
end

delete([tx1 tx2])

% Leave to last - may have to do this explcitily
for inx=1:2
  top(axx(inx))
end

% Last-minute arrangements
delete(yl([1 2 4 5]))
fig2print(gcf,'landscape')
figdisp([],[],[],1)

% To make a decent measurement for visual alignment
% sav=4; dels=[1 2 3 5]
% delete(ha(dels))
% shrink(ha(sav),1/3,1/3)
% % Move to absolute center, for picking
% layout(ha(sav),[1/2 1/2],'mm','xy')
% set(ha(sav),'xtickl',xtox,'xgrid','on')


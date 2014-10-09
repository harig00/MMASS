function mermaid03
% MERMAID03
%
% Simons, Nolet et al., JGR, 2009, FIGURE 5.
%
% Plots the results of positive detections
%
% See also: DETECTS, SIGNALS, WTEVENTS3
% 
% Last modified by fjsimons-at-alum.mit.edu, 1/16/2009

dirs={'11042003','11042003','11042003',...
      '09102004','08092007'};
files={'mpilot1_21118','mpilot1_98947','mpilot1_138106',...
       'mpilot2_79891','mpilot3_4759'};
filts={'lowpass','bandpass','bandpass','bandpass','bandpass'};
fops={{1,2,2},{1,8,2,2},{1,8,2,2},{1,8,2,2},{1,8,2,2}};

% Some arbitrary y-axis scalings
scax=[1e3 1e4 1e4 1e3 5e4];
scay=[1 1 1.2 1 1.5];
% Spectrogram parameters
wsec=10; % This is NOT what I used to do as a default, but it works
wolap=0.875;
% Scalogram parameters
tipo='CDF';
nvm=[2 4];
nd=5;
pph=4;
intel=0;
thresh=1;
div=2.5;
% Now it's the same as in MERMAID13 - soft, div=1, threshold2
div=1;

clf
[ah,ha]=krijetem(subnum(5,4));

for inx=1:5
  disp(' ')
  dir=fullfile('/home/fjsimons/MERMAID',dirs{inx},'EVENTS');
  
  % Time-domain plots of the identified events %%%%%%%%%%%%%%%%%%%%
  axes(ha(inx))
  [x,h,Fs,p(inx),xl(inx),yl(inx),tl(inx),pt01(inx,:),xf]=...
      timdomplot(fullfile(dir,sprintf('%s.sac',files{inx})),2,...
		 filts{inx},fops{inx});
  % Adjust axes to be pretty
  xtix=[-300:150:300];
  % txl=-290;
  % txr=290;
  txl=xtix(1)+range(xtix)/60;
  txr=xtix(end)-range(xtix)/60;;
  set(p(inx),'ydata',get(p(inx),'ydata')/scax(inx))
  set(pt01(inx,1),'ydata',get(pt01(inx,1),'ydata')/scax(inx))
  set(pt01(inx,2),'ydata',get(pt01(inx,2),'ydata')/scax(inx))  
  set(ha(inx),'ylim',get(pt01(inx,2),'ydata')/scay(inx))
  set(ha(inx),'xtick',h.T0+xtix,'xtickl',xtix,'xlim',h.T0+[xtix(1) xtix(end)])
  bottom(pt01(inx,1),ha(inx))
  bottom(pt01(inx,2),ha(inx))
  % Now add magnitude, event depth and distance information on these plots
  ylo=ylim;
  leg1=sprintf(' %s = %3.1f%s','\Delta',h.GCARC,str2mat(176));
  leg2=sprintf('mag = %3.1f ',h.MAG);
  legD=sprintf('d = %i km',round(h.EVDP));
  % And the filter parameters!
  if strcmp(filts{inx},'bandpass')
    leg3=sprintf('%i-%i Hz',fops{inx}{1},fops{inx}{2});
  elseif strcmp(filts{inx},'lowpass')
    leg3=sprintf('%i-%i Hz',0,fops{inx}{1});
  else
    % Highpass? Provide a line here in case necessary
    disp('Do something here')
  end
  
  tx1(inx)=text(h.T0+txr,ylo(1)+range(ylo)/10,leg1,...
		'horizontala','right');
  tx2(inx)=text(h.T0+txr,ylo(2)-range(ylo)/10,leg2,...
		'horizontala','right');
  tx3(inx)=text(h.T0+txl,ylo(2)-range(ylo)/10,leg3,...
		'horizontala','left');
  txD(inx)=text(h.T0+txl+0,ylo(1)+range(ylo)/10,legD,...
		'horizontala','left');
  
  % Should I be adding the STA/LTA results to here as well? Or work with
  % longer series? Or work with hourly sections? or not at all?
  STA=10; LTA=100; TR=2; DTR=1; 
  PEM=100; PET=100; PNL=500; ATL=20;
  % Work on the DATA THAT ARE BEING PLOTTED!! THE PRECISE PARAMETERS OF
  % THE TRIGGERING ALGORITHM WILL HAVE TO BE SET BEFORE EVERY EXPERIMENT
  % AFTER A PHASE OF TESTING... WE HOPE... FOR THE FIRST TWO EXPERIMENTS,
  % WE HAD GOOD RESULTS WITH TWICE DECIMATED, UNFILTERED DATA, BUT WHO
  % KNOWS. IN THE PAPER MAYBE INCLUDE A TABLE
  %   [trigt,stav,ltav,ratio,tim1,tim2,tim3,trigs,dtrigs]=...
  %       stalta(xf,h.DELTA,[h.B h.E],...
  % 	     STA,LTA,TR,DTR,PEM,PET,PNL,ATL);
  %   if ~isempty(trigt)
  %     % This is assuming only two points get triggered
  %     hold on
  %     pst(inx,:)=plot(repmat(trigt,2,1),ylim);
  %     hold off
  %   end
  % In the end, I decided against plotting this in here

  % Spectrogram of the identified events %%%%%%%%%%%%%%%%%%%%
  axes(ha(inx+5))
  % The desired window length, in samples
  wlen=floor(wsec*Fs);
  % The number of frequencies, ideally the length of the window
  nfft=max(2^nextpow2(wlen),512);
  nfft=512; % Actually, keep it the same
  [p2(inx),xl2(inx),yl2(inx),bm{inx},Bl10,F,T]=...
      timspecplot(x,h,nfft,Fs,wlen,wolap,h.B);
  set(ha(inx+5),'xtick',h.T0+xtix,'xtickl',xtix,'xlim',...
		h.T0+[xtix(1) xtix(end)],'ytick',[0:2:10])
  ylo=ylim;
  leg7=sprintf('%i s / %i%s ',wsec,round(100*wolap),'%');
  tx7(inx)=text(h.T0+txr,ylo(2)-range(ylo)/10,leg7,...
		'horizontala','right','FontS',8);
  hold on
  fb7(inx)=fillbox(ext2lrtb(tx7(inx),1,0.7),'w');
  hold off
  set(fb7(inx),'EdgeC','w')
  top(tx7(inx),gca)

  % Scalogram of the identified events %%%%%%%%%%%%%%%%%%%%
  axes(ha(inx+2*5))
  [p3{inx},xl3(inx),yl3(inx),ptt,rgb(inx)]=...
      scalogramplot(x,h,tipo,nvm,nd,pph,intel,thresh,div);
  delete(ptt)
  set(ha(inx+2*5),'xtick',h.T0+xtix,'xtickl',xtix,'xlim',h.T0+[xtix(1) xtix(end)])

  % Now add wavelet information on these plots
  ylo=ylim;  % If not biorthogonal, only provide one slot
  leg4=sprintf('%s(%i,%i) ',tipo,nvm(1),nvm(2));
  if thresh==1
    leg5=sprintf(' %s %3.1f','soft',div);
  end
  tx4(inx)=text(h.T0+txr,ylo(1)+range(ylo)/10,leg4,...
		'horizontala','right','FontS',8);
  tx5(inx)=text(h.T0+txl,ylo(1)+range(ylo)/10,leg5,...
		'horizontala','left','FontS',8);
  % Provide white box behind it - do Fontsize before you do this
  % Ratio to scale white box underneath tex 
  hold on
  fb4(inx)=fillbox(ext2lrtb(tx4(inx)),'w');
  fb5(inx)=fillbox(ext2lrtb(tx5(inx)),'w');
  hold off
  set(fb4(inx),'EdgeC','w')
  top(tx4(inx),gca)
  set(fb5(inx),'EdgeC','w')
  top(tx5(inx),gca)
  % Spectral density of the 1000 s data stream %%%%%%%%%%%%%%%%%%%%
  axes(ha(inx+3*5))
  % lwin=floor(h.NPTS/2);
  lwin=floor(500/h.DELTA);
  olap=70;
  % nfft2=max(2^nextpow2(wlen),512);
  nfft2=512; % Actually, keep it the same
  [p4(inx,:),xl4(inx),yl4(inx)]=specdensplot(x,nfft2,Fs,lwin,olap,1);
  % Make it all relative to zero as in decibels - note that this is
  % relative to ALL of the frequencies and that you later only plot F(2:end)
  wats=get(p4(inx,:),'Ydata'); wats=max(max(cat(1,wats{1:3})));
  % DO NOT DO THIS ANYMORE as per the REVIEWER'S COMMENT
  wats=0;
  for onx=1:4
    set(p4(inx,onx),'Ydata',get(p4(inx,onx),'Ydata')-wats);
  end
  axis tight
  ylo=ylim;
  xlax=xlim;
  % Round off frequency to 10 and redo labels
  Fmax=10;
  set(gca,'xlim',[xlax(1) Fmax])
  mima=[F(2) Fmax];
  poslab=10.^[-3:3];
  poslab=poslab(poslab>=mima(1) & poslab<=mima(2));
  set(gca,'Xtick',poslab,'xtickl',poslab);
  leg6=sprintf('%i s / %i%s',round(lwin/Fs),olap,'%');
  leg6=sprintf('%i s ',round(lwin/Fs));
  tx6(inx)=text(F(2)+(F(3)-F(2))/10,ylo(1)+range(ylo)/10,leg6,...
		'horizontala','left','FontS',8);
  hold on
  pg{inx}=plot(repmat(poslab,2,1),ylo,'k:');
  hold off
end

% Move everything up to make room for the colorbars
movev(ha,.1)
set(p,'Color','k')
%set(pst(~~pst),'Color','w')
%set(pt01,'LineS','-','Color',grey,'LineW',1)
delete(pt01)
for inx=1:4
  serre(ha([1:5]+(inx-1)*5),1/2,'down')
end
moveh(ha(6:10),.005)
moveh(ha(16:end),.005)
% Cosmetic adjustments
% Mark the P-wave arrival time
for inx=1:3
  axx(inx)=xtraxis(ha(1+(inx-1)*5),get(ha(1+(inx-1)*5),'xtick'),...
	      {'' '' 'P' '' ''});
end
[axx(4),xl5,yl5]=xtraxis(ha(16),[10 1 0.1],[10 1 0.1],'period (s)');
for inx=2:5
  axx(3+inx)=xtraxis(ha(15+inx),[10 1 0.1],[]);
end
% Watch out the next line contains a very dangerous hardwire
set(axx(4:end),'Xdir','rev','xlim',[1/Fmax 1/F(2)])

longticks([ah axx])
set([ha xl yl xl2 yl2 xl3 yl3 xl4 yl4 xl5 axx tx1 tx2 tx3 tx4 txD],...
    'FontS',8)

set(p,'LineW',0.25)
set(p4(:,1),'LineW',0.5,'Color','r');
set(p4(:,[2 3]),'LineW',0.5,'Color',grey);
set(p4(:,4),'MarkerS',2,'Marker','o','MarkerF','r','MarkerE','r');
delete(xl(1:4))
delete(xl2(1:4))
delete(xl3(1:4))
delete(xl4(1:4))
set(xl(end),'string','time (s)')
set(yl,'string','sound pressure (scaled hydrophone counts)')
set(xl2(end),'string','time (s)')
set(xl3(end),'string','time (s)')
set(xl4(end),'string','frequency (Hz)')
set(yl2,'string','frequency (Hz)')
set(yl3,'string','scale')
set(yl4,'string','log spectral density')
for inx=1:4
  nolabels(ha([1:4]+(inx-1)*5),1) % This must be first
end
delete(tl)
[bh,th]=label(ha(1:10),'ll',10,[],[],[],[],0.75);
[bh(11:15),th(11:15)]=label(ha(11:15),'ul',10,10,[],[],[],0.75);

% Now get the extent of the added labels and shift legD by that amount
for inx=1:5
  shifit=diff(unique(get(bh(inx),'XData')));
  set(txD(inx),'Position',getpos(txD(inx))+[1.4*shifit 0 0])
end

% Get rid of them
delete([bh(6:end) th(6:end)])

colormap jet 
axes(ha(10))
cb(1)=colorbarf('hor',8,'Helvetica',...
	[getpos(ha(10),[1 2 3]) 0.015]+[0 -0.07 0 0]);

axes(ha(15))
cb(2)=colorbarf('hor',8,'Helvetica',...
	[getpos(ha(15),[1 2 3]) 0.015]+[0 -0.07 0 0]);

% Very important adjustment to the color axis
for inx=1:5
  cdt=get(p2(inx),'CData');
  k=mean(cdt(:));
  l=std(cdt(:));
  levs=3;
  set(ha(5+inx),'Clim',[k-levs*l k+levs*l])
end
xcb1=get(cb(1),'xlim');
axes(cb(1))
xcb(1)=xlabel('log spectral density (stdev from mean)','FontS',8);
set(cb(1),'xtick',linspace(xcb1(1),xcb1(2),2*levs+1),...
	  'xtickl',[-levs:1:levs])

% Very important adjustment to the color axis assuming indeed this goes
% to the hardwired three times the standard deviation!!
axes(cb(2))
set(gcf,'NextPlot','add')
image(repmat([10:-0.5:0]/10,[1 1 3]))
xcb2=get(cb(2),'xlim');
set(cb(2),'xtick',linspace(xcb2(1),xcb2(2),4),...
	  'xtickl',[0 1 2 3],'ytick',[])
xcb(2)=xlabel('wavelet coeff magnitude (stdev)','FontS',8);
longticks(cb,2)

% The wavelet legend is everywhere the same - cut
delete([tx4(:) ;  fb4(:)])
% And so are these:
delete([tx7(:) ;  fb7(:)])
delete([tx5(:) ;  fb5(:)])
delete([tx6(:)])

% Leave to last - may have to do this explcitily
for inx=1:4
  top(axx(inx))
end

% Last-minute arrangements
delete(yl([1 2 4 5]))
fig2print(gcf,'landscape')
figdisp([],[],[],0)
%disp('!degs /home/fjsimons/EPS/mermaid03.eps')

function mermaid05
% MERMAID05
%
% Simons, Nolet et al., JGR, 2009, FIGURE 6.
%
% Plots the best estimate of the noise spectrum of the three experiments
% (at depth, no events nor pumps)... using SPECDENSPLOT and PCHAVE.
%
% Last modified by fjsimons-at-alum.mit.edu, 01/16/2009

% Specify the experiment dat
dirs={'11042003','09102004','08092007'};

% Specify within those directories the running number of the hourlong
% data sections with the 10-minute overlap. Do this from
% file:///home/fjsimons/PIX/MERMAID/11042003/HRSECTIONS/web/index.html
% file:///home/fjsimons/PIX/MERMAID/09102004/HRSECTIONS/web/index.html
% file:///home/fjsimons/PIX/MERMAID/08092007/HRSECTIONS/web/index.html
% and cross-reference with the "meaningful selections"
% file:///home/fjsimons/PIX/MERMAID/11042003/SELECTIONS/web/6.html
% file:///home/fjsimons/PIX/MERMAID/09102004/SELECTIONS/web/6.html
% For the third experiment, I still have to make this "SELECTIONs"
expo{1}=[9 11:13 15:23 28 29 31 32 35:45]; 
expo{2}=[10:34 36 37 39:46];
expo{3}=[9:19 21:25 27 29:30 32 34 36:51 53];

% Prepare the figure panels
[ah,ha]=krijetem(subnum(1,3));

% Go through and do the analysis
for index=1:3
  axes(ah(index))
  for ondex=1:length(expo{index})
    % Read in the data, hourlong sections, with 10 minutes of overlap
    filenam1=fullfile('/home/fjsimons/MERMAID',dirs{index},'HRSECTIONS',...
		      'REALDATA',sprintf('mpilot%i_%3.3i.sac',...
				  index,expo{index}(ondex)));
    % display(sprintf('Reading %s',filenam1))
    [x,h,tl]=readsac(filenam1,0,'l');
    Fs(ondex)=1/h.DELTA;
    npts(ondex)=h.NPTS;
    % display(sprintf('Sampling frequency %5.3e',Fs))
    bigex{ondex}=x;
  end
  % Check the sampling rate and sample size are consistent
  difer(sum(diff(Fs))); Fs=Fs(1);
  difer(sum(diff(npts))); npts=npts(1); 

  % Now do the spectral analysis and make a plot 
  nfft2=512;
  % Single window per segment
  lwin=npts;
  % Any amount of overlap will do in this case
  olap=0;
  [p4(index,:),xl4(index),yl4(index),F]=...
      specdensplot(bigex,nfft2,Fs,lwin,olap,1);
  % Make it all relative to zero as in decibels - note that this is
  % relative to ALL of the frequencies and that you later only plot F(2:end)
  wats=get(p4(index,:),'Ydata'); wats=max(max(cat(1,wats{1:3})));
  % DO NOT DO THIS ANYMORE as per the REVIEWER'S COMMENT
  wats=0;
  for onx=1:4
    set(p4(index,onx),'Ydata',get(p4(index,onx),'Ydata')-wats);
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
  hold on
  pg{index}=plot(repmat(poslab,2,1),ylo,'k:');
  hold off
end

% Cosmetics
shrink(ah,1,4)
for index=1:3
  [axx(index),xl5(index),yl5(index)]=...
      xtraxis(ha(index),[10 1 0.1],[10 1 0.1],'period (s)');
end
set(axx,'Xdir','rev','xlim',[1/Fmax 1/F(2)])
longticks([ah axx])
set([ha xl4 yl4 xl5 axx],'FontS',8)
set(p4(:,1),'LineW',0.5,'Color','r');
set(p4(:,[2 3]),'LineW',0.5,'Color',grey);
set(p4(:,4),'MarkerS',2,'Marker','o','MarkerF','r','MarkerE','r');
set(xl4,'string','frequency (Hz)')
set(yl4,'string','log spectral density')
delete(yl4(2:end))
serre(ah,1/5,'across')
serre(axx,1/5,'across')
% Check out them tricks for the logarithmic labels
[bh,th]=label(ha,'ll',10,[],[],[],1/30,1,1);
moveh(bh,-0.205)
moveh(th,-0.2065)

moveh(ah(1),-.1/3)
moveh(axx(1),-.1/3)
moveh(ah(2),-.1/6)
moveh(axx(2),-.1/6)

for inx=1:3
  top(axx(inx))
end

fig2print(gcf,'portrait')
figdisp([],[],[],0)

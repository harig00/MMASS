function dsmcouplinglog(as,np)
% DSMCOUPLINGLOG(as,np)
%
% Makes a LOG plot of the coupling matrix for boxcar windows
% Dahlen & Simons (2007), Figure 7
%
% INPUT:
%
% as      1 also plot asymptotic relation [default: 0]
% np      Number of panels [default: 4]
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2007

defval('Lmax',150)
defval('as',0)
defval('np',4)

% Bandwidths of left and right column of panels
L1=[5 10];
L2=[20 30];

% X-axis limits
xma1=20;
xma2=40;
xti1=[20 20];
xti2=[40 40];
xtag=['offset from target degree l'' - l'];

% Y-axis ticks and limits 
ytag=sprintf(['M_{ll''}  (dB)']);

clf
nr=3;
[ah,ha]=krijetem(subnum(nr,2));
fig2print(gcf,'portrait')

% FIRST COLUMN
for ondex=1:length(L1)
  axes(ha(ondex))
  
  L=L1(ondex);

  % Get the eigenvalue-weighted multitaper coupling kernel
  K=mcouplings(L,Lmax);
 
  if as==1
    % And get the asymptotic representation fake for the last one
    [Kas,elas]=universal(L,0);
    % What seems to be the scaling factor here?
    scK=K(llast,llast)/max(Kas);
    disp(sprintf('Scaling is %8.3f',scK))
    Kas=Kas*100;
  end
      
  % Somehow normalize this to the effective bandwidth as percentage
  % leaked to, from and within the effective bandwidth
 
  % Which one are we plotting? The last complete degree (not index)
  llast=Lmax-L;
  % Check row sum
  Ksum=sum(100*K(llast+1,max(1,llast-L):Lmax+1));
  disp(sprintf('Row sum at degree shown %8.3f',Ksum))
  difer(Ksum-100)

  % Make a log scale version
  K(K==0)=NaN;
  K=decibel(K(llast+1,:));

  minb=min(K);
  % Actually, adjust to a common minimum scale
  minb=-4;
  maxb=-minb;
  tix=[0:1:maxb];
  widr=[0 0.75];

  if as==1
    plot(elas,Kas*scK,'k','linew',1)
    hold on
  end
  
  b=bar(0-llast:Lmax-llast,K-minb,1);
  set(b,'FaceC',grey,'EdgeC','k')
  xlim([-xma1 xma1])
  
  set(gca,'ytick',sort(maxb-tix),'ytickl',fliplr(-tix))
  ylim([0 max(K-minb)]+widr)
  
  tixx{ondex}=[0:xti1(ondex):xma1 -L L];
  xl(ondex)=xlabel(xtag);
  yl(ondex)=ylabel(ytag);
  % Put in bandwidth labels
  legsi{ondex}=sprintf(' L = %i ',L);
  [bh(ondex),th(ondex)]=boxtex('ur',ha(ondex),legsi{ondex},12,1,1,1,...
			       1.01,1);
end

% SECOND COLUMN
for ondex=1:length(L2)
  panid=ondex+length(L1)+(nr-length(L1));
  axes(ha(panid));
  
  L=L2(ondex);

  % Get the boxcar coupling kernel
  K=mcouplings(L,Lmax);

  if as==1
    % And get the asymptotic representation fake for the last one
    [Kas,elas]=universal(L,0);

    % What seems to be the scaling factor here?
    scK=K(llast,llast)/max(Kas);
    disp(sprintf('Scaling is %8.3f',scK))
    Kas=Kas*100;
  end

  % Somehow normalize this to the effective bandwidth as percentage
  % leaked to, from and within the effective bandwidth

  % Which one are we plotting? 
  llast=Lmax-L;
  % Check the row sum
  Ksum=sum(100*K(llast+1,max(1,llast-L):Lmax+1));
  disp(sprintf('Row sum at degree shown %8.3f',Ksum))
  difer(Ksum-100)
  
  % Make a log scale version
  K(K==0)=NaN;
  K=decibel(K(llast+1,:));
  minb=min(K);
  % Actually, adjust to a common minimum scale
  minb=-6;
  maxb=-minb;
  tix=[0:2:maxb];
  widr=[0 0.75];

  if as==1
    plot(elas,Kas,'k','linew',1)
    hold on
  end
  b=bar(0-llast:Lmax-llast,K-minb,1);
  set(b,'FaceC',grey,'EdgeC','k')
  xlim([-xma2 xma2])

  set(gca,'ytick',sort(maxb-tix),'ytickl',fliplr(-tix))
  ylim([0 max(K-minb)]+widr)
  
  tixx{panid}=[0:xti2(ondex):xma2 -L L];
  xl(panid)=xlabel(xtag);
  yl(panid)=ylabel(ytag);
  % Put in bandwidth labels
  legsi{panid}=sprintf(' L = %i ',L);
  [bh(panid),th(panid)]=...
      boxtex('ur',ha(panid),legsi{panid},12,1,1,...
	     1.01,1);
end

yma1=[15 9];
yma2=[6 3];

% Cosmetic arrangements
longticks(ah)
movs=1e9; % Don't do this here unless you do it everywhere
for ondex=1:length(L1)
  set(ha(ondex),'xtick',unique([-tixx{ondex} tixx{ondex}]),'xgrid','on')
  movev(bh(ondex),-yma1(ondex)/movs)
  movev(th(ondex),-yma1(ondex)/movs)
end
for ondex=1:length(L2)
  panid=ondex+length(L1)+(nr-length(L1));
  set(ha(panid),'xtick',unique([-tixx{panid} tixx{panid}]),'xgrid','on')
  movev(bh(panid),-yma2(ondex)/movs)
  movev(th(panid),-yma2(ondex)/movs)
end
% Four panels or six?
if np==4
  delete(xl([1 4]))
  delete(ah(5:6))
else
  delete(xl([1:2 4:5]))
end

serre(ha(1:2),1/30,'down')
serre(ha(4:5),1/30,'down')

set(gcf,'color','w','inverthardcopy','off')
if as==1
  figdisp([],as)
else
  figdisp
end


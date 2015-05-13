function dsbcoupling(as,offs)
% DSBCOUPLING(as,offs)
%
% Makes a plot of the coupling matrix for boxcar windows
% Dahlen & Simons (2007), Figure 5.
%
% INPUT:
%
% as      1 also plot asymptotic relation [default: 0]
% offs    An offset in the last degree plotted [default: 0]
%
% Last modified by fjsimons-at-alum.mit.edu, 12/22/2006

% Run it twice, to be safe

defval('Lmax',150)
defval('as',0)
defval('offs',0);

TH1=[10 20 30]; 
TH2=[15 10 5];
TH2=[30 20 10];

% EL-spacing
els=2;
els2=1;

% X-tick interval 
tb=[20 10 5];
tb2=[5 2 2];
xtag=['offset from target degree l'' - l'];
% Y-axes
tbb=[10 15 25];
tbb=[5 10 15];
tbb2=[70 80 100];

clf
[ah,ha]=krijetem(subnum(3,2));
fig2print(gcf,'portrait')

% Logo scale and location
logsk=20;
cloc='ll';

aN=3;
for ondex=1:length(TH1)
  axes(ha(ondex))
  % Get effective bandwidth thingies
  for N=1:max(5,aN)
    skel(N)=max(roots([1 1 -(N*180./TH1(ondex))^2]));
  end
    
  % Get the boxcar coupling kernel
  K=bcoupling(TH1(ondex),Lmax,1);

  % Last l-prime we'll be showing today (a degree, not an index)
  llast=Lmax-ceil(skel(aN))-offs;
  disp(sprintf('We show target degree l = %i',llast))
  disp(sprintf('Leftmost portion l = %i',llast-ceil(skel(aN))))

  % Check the row sum, which isn't going to be 100, but just about
  Ksum=sum(K(llast+1,...
	     max(1,llast-ceil(skel(aN)))+1:min(llast+ceil(skel(aN))+1,Lmax+1)));
  disp(sprintf('Row sum at degree shown %8.3f',Ksum*100))

  if as==1
    % And get the asymptotic representation fake for the last one
    [Kas,elas]=buniversal(TH1(ondex),Lmax,0,1);

    % What seems to be the scaling factor here?
    scK=K(llast,llast)/max(Kas);
    disp(sprintf('Scaling is %8.3f',scK))
    Kas=Kas*100;
    plot(elas,Kas*scK,'k','linew',1)
    hold on
  end
  % Somehow normalize this to the effective bandwidth as percentage
  % leaked to, from and within the effective bandwidth
  K=K*100;

  b=bar(0-llast:Lmax-llast,K(llast+1,:),1);
  set(b,'FaceC',grey,'EdgeC','k')
  ylim([0 tbb(ondex)])
  xlim([-skel(aN) skel(aN)])
  tix{ondex}=0:tb(ondex):round(skel(aN)/10)*10;
  xl(ondex)=xlabel(xtag);
  yl(ondex)=ylabel(sprintf(['100 %s K_{ll''}  (%s)'],'\times','%'));
  % Put in cap size labels
  legsi{ondex}=sprintf('  %s = %i%s','\Theta',TH1(ondex),str2mat(176));
  [bh(ondex),th(ondex)]=boxtex('ur',ha(ondex),legsi{ondex},12,1,1);
  [bh2(ondex),th2(ondex)]=boxtex('ul',ha(ondex),...
			 sprintf('%i%s',100-round(Ksum*100),'%'),12,1,1);
  delete(bh2(ondex))
  % Equivalent wavelength axis
  [ax(ondex),axl(ondex)]=...
      xtraxis(ha(ondex),sort([-skel(1:aN) 0 skel(1:aN)]),-aN:aN,...
	      ['number of wavelengths per cap']); 
  set(ax(ondex),'xgrid','on')

  % Put the cap logo on
  lah(ondex)=caplogo(ha(ondex),1,cloc,1/logsk,1/logsk,20);
end

clear llast
% Now for the double cap
for ondex=1:length(TH2)
  axes(ha(ondex+length(TH1)))
  % Get effective bandwidth thingies
  for N=1:max(5,aN)
    skel(N)=max(roots([1 1 -(N*180./(90-TH2(ondex)))^2]));
  end
  
  % Get the boxcar coupling kernel
  K=bcoupling(TH2(ondex),Lmax,2);

  % Target degree (not an index)
  llast=Lmax-ceil(skel(aN))-offs;
  disp(sprintf('We show target degree l = %i',llast))

  % Check the row sum, which isn't going to be 100, but just about
  Ksum=sum(K(llast+1,...
	     max(1,llast-ceil(skel(aN)))+1:min(llast+ceil(skel(aN))+1,Lmax+1)));
  disp(sprintf('Row sum at degree shown %8.3f',Ksum*100))

  if as==1
    % And get the asymptotic representation fake for the last one
    [Kas,elas]=buniversal(TH2(ondex),Lmax,0,2);
    
    % What seems to be the scaling factor here?
    scK=K(llast,llast)/max(Kas);
    disp(sprintf('Scaling is %8.3f',scK))
    Kas=Kas*100;
    % Plot asymptotic result?
    plot(elas,Kas,'k','linew',1)
    hold on
  end

  % Somehow normalize this to the effective bandwidth as percentage
  % leaked to, from and within the effective bandwidth
  K=K*100;

  b=bar(0-llast:Lmax-llast,K(llast+1,:),1);
  set(b,'FaceC',grey,'EdgeC','k')
  
  ylim([0 tbb2(ondex)])
  xlim([-skel(aN) skel(aN)])
  tix{ondex+length(TH1)}=0:tb2(ondex):round(skel(aN)/10)*10;
  xl(ondex+length(TH1))=xlabel(xtag);
  yl(ondex+length(TH1))=ylabel(sprintf(['100 %s K_{ll''}  (%s)'],'\times','%'));
  % Put in cap size labels
  legsi{ondex+length(TH1)}=...
      sprintf('  %s = %i%s','\Theta',90-TH2(ondex),str2mat(176));
  [bh(ondex+length(TH1)),th(ondex+length(TH1))]=...
      boxtex('ur',ha(ondex+length(TH1)),legsi{ondex+length(TH1)},12,1,1);
  [bh2(ondex+length(TH1)),th2(ondex+length(TH1))]=...
      boxtex('ul',ha(ondex+length(TH1)),...
	     sprintf('%i%s',100-round(Ksum*100),'%'),12,1,1);
  delete(bh2(ondex+length(TH1)))

  % Equivalent wavelength axis
  [ax(ondex+length(TH1)),axl(ondex+length(TH1))]=...
      xtraxis(ha(ondex+length(TH1)),sort([-skel(1:aN) 0 skel(1:aN)]),...
	      -aN:aN,['number of wavelengths per cap']); 
  set(ax(ondex+length(TH1)),'xgrid','on')

  % Put the cap logo on
  lah(ondex+length(TH1))=...
      caplogo(ha(ondex+length(TH1)),3,cloc,1/logsk,1/logsk,20);
end

% Cosmetic arrangements
longticks([ah ax])
delete(axl([2:3 5:6]))
delete(xl([1:2 4:5]))
nolabels(ax([2:3 5:6]),1)
for ondex=1:length(ha)
  set(ha(ondex),'xtick',unique([-tix{ondex} tix{ondex}]))
end
set(th2,'FontS',8)
movev(lah,.003)

set(gcf,'color','w','inverthardcopy','off')
if as==1
  figdisp([],as)
else
  figdisp
end

% Make the extra axis come last
for index=1:length(ah)
  axes(ah(index))
  set(ah(index),'color','none')
end

function dsbcouplings
% DSBCOUPLINGs
%
% Makes a cascade plot of the coupling matrix for boxcar windows
% Dahlen & Simons (2007), Figure 4
%
% Last modified by fjsimons-at-alum.mit.edu, 12/21/2006

defval('Lmax',100)

TH1=[10 20 30]; 
TH2=[30 20 10];

% Number, not spacing, of degrees 
els=[3 4 4];
els2=[3 4 4];

% X-tick interval 
tb=20;
tb2=20;
% Y-axes
fax=1.3;
tbb=fax*[10 17.5 25];
tbb2=fax*[60 80 100];
% Y-tick interval
ty=[3 5 10];
tyt=[9 15 20];
ty2=[20 40 50];
tyt2=[60 80 100];

clf
[ah,ha]=krijetem(subnum(3,2));

% Logo scale and location
logsk=20;
cloc='ll';

aN=3;
% This for the SINGLE CAP
for ondex=1:length(TH1)
  axes(ha(ondex))
  % Get effective bandwidth thingies
  for N=1:5
    skel(N)=max(roots([1 1 -(N*180./TH1(ondex))^2]));
  end

  K=bcoupling(TH1(ondex),Lmax,1);
  % Somehow normalize this to the effective bandwidth as percentage
  % leaked to, from and within the effective bandwidth
  K=K*100;
 
  clear b
  % Do it until Lmax minus the point at which the power spectrum
  % essentially drops to zero so you don't run into truncation effects
  numax=0;
  ax(1,ondex)=ha(ondex);
  % Make smaller and move down in anticipation 
  shrink(ax(1,ondex),1.35,1.35)
  movev(ax(1,ondex),-0.05)
  % for index=floor(linspace(0,Lmax-skel(aN),els(ondex)))
  for index=cumsum([0 repmat(20,1,els(ondex)-1)])
    if numax>0
      ax(numax+1,ondex)=laxis(ax(numax,ondex),...
			      -(1/6)*[min(els)-1]/[els(ondex)-1],...
			      1/5);
    end
    numax=numax+1;
    b(numax)=bar(0:Lmax,K(index+1,:),1);
    set(b(numax),'FaceC',grey)
    [i,j]=max(K(index+1,:));
    t(numax,ondex)=text(j-1,i+tbb(ondex)/15,sprintf('%i',round(i)),...
			'horizontala','center');
    t2(numax,ondex)=text(102,tbb(ondex)/10,sprintf('l = %i',round(index)),...
			'horizontala','right');
    if numax==1
      xl(ondex)=xlabel('degree l''');
      yl(ondex)=ylabel(sprintf(['100 %s K_{ll''}  (%s)'],'\times','%'));
    end
  end
% set(b(end),'FaceC','w')
  % Reorder axes
  for index=numax:-1:1
    axes(ax(index,ondex))
    set(ax(index,ondex),'ylim',[0 tbb(ondex)])
    set(ax(index,ondex),'ytick',[0:5*(1+[ondex==3]):tbb(ondex)/fax])
    set(ax(index,ondex),'xlim',[-2 100])
  end
  set(ax(2:end,ondex),'ycolor','w')
  set(ax(1:end,ondex),'color','none')
  nolabels(ax(2:end,ondex),1)
  % Put the cap logo on
  lah(ondex)=caplogo(ha(ondex),1,cloc,1/logsk,1/logsk,20);
end
moveh(lah,0.29)
movev(lah,-0.025)
for ondex=1:3
  % Put in cap size labels
  legsi{ondex}=...
      sprintf('  %s = %i%s','\Theta',TH1(ondex),str2mat(176));
  [bh(ondex),th(ondex)]=...
      boxtex('ul',ha(ondex),legsi{ondex},10,[],0.7,0.95);
end

set(ax(~~ax),'box','off')
set(ax(~~ax),'xtick',[0:tb:100])
for index=1:length(TH1)
  set(ax(1,index),'ytick',[0:ty(index):tyt(index)])
end
longticks(ax(~~ax))
set(t(~~t),'FontS',7)
set(t2(~~t2),'FontS',9)
set(xl(~~xl),'FontS',9)
set(yl(~~yl),'FontS',9)

clear ax
% Now for the double cap
for ondex=1:length(TH2)
  axes(ha(ondex+length(TH1)))
  % Get effective bandwidth thingies
  for N=1:5
    skel(N)=max(roots([1 1 -(N*180./(90-TH2(ondex)))^2]));
  end
  
  % Get coupling kernel
  K=bcoupling(TH2(ondex),Lmax,2);
  
  % Somehow normalize this to the effective bandwidth as percentage
  % leaked to, from and within the effective bandwidth
  K=K*100;
 
  clear b 
  % Do it until Lmax minus the point at which the power spectruym
  % essentially drops to zero so you don't run into truncation effects
  numax=0;
  ax(1,ondex)=ha(ondex+length(TH1));
  % Make smaller and move down in anticipation 
  shrink(ax(1,ondex),1.35,1.35)
  movev(ax(1,ondex),-0.05)
  % for index=floor(linspace(0,Lmax-skel(aN),els(ondex)))
  for index=cumsum([0 repmat(20,1,els2(ondex)-1)])
    if numax>0
      ax(numax+1,ondex)=laxis(ax(numax,ondex),...
			      -(1/6)*[min(els2)-1]/[els2(ondex)-1],...
			      1/5);
    end
    numax=numax+1;
    b(numax)=bar(0:Lmax,K(index+1,:),1);
    set(b(numax),'FaceC',grey)
    [i,j]=max(K(index+1,:));
    t(numax,ondex)=text(j-1,i+tbb2(ondex)/15,sprintf('%i',round(i)),...
			'horizontala','center');
    t2(numax,ondex)=text(102,tbb2(ondex)/10,sprintf('l = %i',round(index)),...
			'horizontala','right');
    if numax==1
      xl(ondex)=xlabel('degree l''');
      yl(ondex)=ylabel(sprintf(['100 %s K_{ll''}  (%s)'],'\times','%'));
    end
  end
% set(b(end),'FaceC','w')
  % Reorder axes
  for index=numax:-1:1
    axes(ax(index,ondex))
    set(ax(index,ondex),'ylim',[0 tbb2(ondex)])
    set(ax(index,ondex),'xlim',[-2 100])
  end
  set(ax(2:end,ondex),'ycolor','w')
  set(ax(1:end,ondex),'color','none')
  nolabels(ax(2:end,ondex),1)
  % Put the cap logo on
  lah(ondex+length(TH1))=...
      caplogo(ha(ondex+length(TH1)),3,cloc,1/logsk,1/logsk,20);
end
moveh(lah(4:6),0.29)
movev(lah(4:6),-0.025)
set(t(1,1),'horizontala','left')
set(t(1:2,2:3),'horizontala','left')

for ondex=1:3
  % Put in cap size labels
  legsi{ondex+length(TH1)}=...
      sprintf('  %s = %i%s','\Theta',90-TH2(ondex),str2mat(176));
  [bh(ondex+length(TH1)),th(ondex+length(TH1))]=...
      boxtex('ul',ha(ondex+length(TH1)),legsi{ondex+length(TH1)},10,[],0.7,0.95);
end

set(ax(~~ax),'box','off')
set(ax(~~ax),'xtick',[0:tb2:100])
for index=1:length(TH2)
  set(ax(1,index),'ytick',[0:ty2(index):tyt2(index)])
end
longticks(ax(~~ax))
set(t(~~t),'FontS',7)
set(t2(~~t2),'FontS',9)
set(xl(~~xl),'FontS',9)
set(yl(~~yl),'FontS',9)

set(ah,'fonts',9)

set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')
figdisp

moveh(bh,-1)
for ondex=1:3
  movev(bh(ondex),tbb(ondex)/fax/40)
  movev(th(ondex),tbb(ondex)/fax/40)
end
for ondex=1:3
  movev(bh(ondex+length(TH1)),tbb2(ondex)/fax/40)
  movev(th(ondex+length(TH1)),tbb2(ondex)/fax/40)
end

disp('Run again in same figure window to get everything right!')

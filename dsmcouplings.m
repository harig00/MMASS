function dsmcouplings
% DSMCOUPLINGS
%
% Cascade of coupling matrix for eigenvalue-weighted multitaper windows 
% Dahlen & Simons (2007), Figure 6
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2007

defval('Lmax',150)

% Bandwidths of left and right column of panels
L1=[10 20];

% Number, not spacing, of degrees 
els=[6 6];
elsp=10;

% X-tick interval 
tb=20;
xma=100;
% Y-axes
fax=1.3;
tbb=fax*[19.5 10.5];
xtag=['degree l'''];
% Y-tick interval
ty=[5 3];
tyt=[20 9];
ytag=sprintf(['100 %s M_{ll''}  (%s)'],'\times','%');

clf
[ah,ha]=krijetem(subnum(3,2));

% This for the FIRST PANEL
for ondex=1:length(L1)
  axes(ha(ondex))

  % Get coupling kernel
  K=mcouplings(L1(ondex),Lmax);
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
  movev(ax(1,ondex),-0.1)
  for index=cumsum([0 repmat(elsp,1,els(ondex)-1)])
    if numax>0
      ax(numax+1,ondex)=laxis(ax(numax,ondex),...
			      -(1/3)*[min(els)-1]/[els(ondex)-1],...
			      1/5);
    end
    numax=numax+1;
    b(numax)=bar(0:Lmax,K(index+1,:),1);
    set(b(numax),'FaceC',grey,'EdgeC','k')
    [i,j]=max(K(index+1,:));
    % Put maximum value right on top of the peak of the histogram
    t(numax,ondex)=text(j-1,i+tbb(ondex)/15,sprintf('%i',round(i)),...
			'horizontala','center');
    % Put the degree to the right of the individual axes
    t2(numax,ondex)=text(xma+xma/50,tbb(ondex)/10,...
			 sprintf('l = %i',round(index)),...
			'horizontala','right');
    if numax==1
      xl(ondex)=xlabel(xtag);
      yl(ondex)=ylabel(ytag);
    end
  end

  % Reorder axes
  for index=numax:-1:1
    axes(ax(index,ondex))
    set(ax(index,ondex),'ylim',[0 tbb(ondex)])
    set(ax(index,ondex),'ytick',[0:5*(1+[ondex==3]):tbb(ondex)/fax])
    set(ax(index,ondex),'xlim',[-2 xma])
  end
  set(ax(2:end,ondex),'ycolor','w')
  set(ax(1:end,ondex),'color','none')

  nolabels(ax(2:end,ondex),1)
end

for ondex=1:length(L1)
  % Put in bandwidth labels
  legsi{ondex}=sprintf('L = %i ',L1(ondex));
  [bh(ondex),th(ondex)]=...
      boxtex('ul',ha(ondex),legsi{ondex},10,[],0.7,0.95);
end

set(ax(~~ax),'box','off')
set(ax(~~ax),'xtick',[0:tb:xma])
for index=1:length(L1)
  set(ax(1,index),'ytick',[0:ty(index):tyt(index)])
end
longticks(ax(~~ax))
set(t(~~t),'FontS',7)
set(t2(~~t2),'FontS',9)
set(xl(~~xl),'FontS',9)
set(yl(~~yl),'FontS',9)

longticks(ax(~~ax))
set(t(~~t),'FontS',7)
set(t2(~~t2),'FontS',9)
set(xl(~~xl),'FontS',9)
set(yl(~~yl),'FontS',9)

set(ah,'fonts',9)

set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')

moveh(bh,-1)
for ondex=1:length(L1)
  movev(bh(ondex),tbb(ondex)/fax/40)
  movev(th(ondex),tbb(ondex)/fax/40)
end

delete(ha(3:6))

figdisp

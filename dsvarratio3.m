function dsvarratio3
% DSVARRATIO3
%
% Plots the variance ratio of the eigenvalue-weighted multitaper variance
% to the whole-sky estimate
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2007

% LATER ON REWRITE THIS TO USE THE NEW FUNCTION MVARRATIOS

Lmax=50;
l=0:Lmax;

% Various bandwidths to be plotted down the panel rows
EL=[10 20];

% Single-cap parameters per bandwidth
TH{1}=[15 20 30 60];
TH{2}=[15 20 30 60];

% Double-cap parameters per bandwidth, stuff to subtract
TH2{1}=90-[40 50 60 80];
TH2{2}=90-[40 50 60 80];

% The maximum value on the y-axis
ymax=[1.8 1.4];
ymax2=[0.22 0.13];

clf
% Create figure panels
[ah,ha]=krijetem(subnum(length(EL),2));

% Must do this before the plot since the caplogo comes last
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')
serre(ah(1:2),1/4,'across')
serre(ah(3:4),1/4,'across')

% Make the panels
[xl(1),yl(1),px{1},pl(1),p{1}]=...
    doit(EL(1), TH{1},1,Lmax,ah(1),ymax(1),l);
[xl(2),yl(2),px{2},pl(2),p{2}]=...
    doit(EL(1),TH2{1},2,Lmax,ah(2),ymax2(1),l);
[xl(3),yl(3),px{3},pl(3),p{3}]=...
    doit(EL(2), TH{2},1,Lmax,ah(3),ymax(2),l);
[xl(4),yl(4),px{4},pl(4),p{4}]=...
    doit(EL(2),TH2{2},2,Lmax,ah(4),ymax2(2),l);

% And plot
set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')

% What's shown here is the empirical rule with the 0.88 exponent which
% works at all degrees, as soon as A is big enough
%delete(pl)
delete(xl([1 2]))
delete(yl([2 4]))

nolabels(ah(1:2),1)

figdisp([],[],[],0)
!degs /home/fjsimons/EPS/dsvarratio3.eps

disp('Rerun if it doesn''t look quite right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xl,yl,px,pl,p]=doit(L,TH,sord,Lmax,ah,ymax,l)
% Note it ALWAYS does the essential checks even with this to 0
defval('xver',0)
% Get all the ZEROJ coefficients at the same time... not a huge timesaver
[allW,C0,S0,Leff]=zeroj(repmat(0:2:2*L,1,Lmax+1),...
	   gamini(0:Lmax,L+1),gamini(0:Lmax,L+1));

% The next verification can be slow and is rarely necessary
if xver==1 & 1==3
  difer(allW-threej(repmat(0:2:2*L,1,Lmax+1),...
	   gamini(0:Lmax,L+1),gamini(0:Lmax,L+1)));
end

% Calculate the matrix that goes into this - this takes most of the time
% Always only get the evens since we're studying l=l, the variance
[Gp,p,K]=gammap(L,TH,sord,1,1);

% Make the whole-sphere (WS) approximation: NOTE THIS IS THE A=4pi
% MULTITAPER and not the UNTAPERED WS result
bigS=gamini([0:L],(L+1))';
bigSp=repmat([0:L]',(L+1),1);
GpWS=repmat(NaN,1,2*L);
% Rather, keep the ones that you had already!
if Leff<2*L
  [jk,C0,S0]=zeroj(0,0,0,2*L);
  Leff=2*L;
end

for pWS=0:2:2*L
  GpWS(pWS+1)=sum((2*bigS'+1).*(2*bigSp'+1).*...
		 zeroj(bigS,pWS,bigSp,Leff,[],C0,S0).^2);
end
% This must be exact when compared to gammap(4*pi)
GpWS=GpWS(1:2:end)*4*pi/(L+1)^4;
if xver==1
  % Under the single cap thing: single cap of entire globe
  difer(GpWS-gammap(L,180,1,1,1));
  disp('Check for WHOLE-SPHERE from single cap passed')
  % Under the double cap thing: subtract belt of nothing
  difer(GpWS-gammap(L,0,2,1,1));
  disp('Check for WHOLE-SPHERE from double cap passed')
end

% Better focus exclusively on those for which K>1 at least
Gp=Gp(K>=1,:); TH=TH(K>=1); K=K(K>=1);

% Only now initialize v and its approximations
v=repmat(NaN,length(TH),length(l));
vll=repmat(NaN,length(TH),1);
vWS=repmat(NaN,1,length(l));

% Better get all of the wigner0j symbols at once here
for ixl=1:length(l)
  % Now we're here, we can do slightly more right away:
  if xver==1
    % Calculate and verify
    [W,pp]=wigner0j(2*L,l(ixl),l(ixl));
    difer(p-pp(1:2:end))
    difer(allW((L+1)*l(ixl)+1:(L+1)*(l(ixl)+1))-W(1:2:end))
    W=W(1:2:end);
  else
    % Stick with the one-blow precalculated ones
    W=allW((L+1)*l(ixl)+1:(L+1)*(l(ixl)+1));
  end
  % Only select the evens since we're doing variance at equal l=l'
  % And calculate the multitaper covariance ratio
  % For the single cap
  v(:,ixl)=(2*l(ixl)+1)/(4*pi)*[repmat(2*p+1,length(TH),1).*Gp]*[W.^2]'; 
  vWS(ixl)=(2*l(ixl)+1)/(4*pi)*[repmat(2*p+1,1,1).*GpWS]*[W.^2]'; 
end

if xver==1
  % What should the zero-l limit be?
  bigS=gamini([0:L],(L+1))';
  bigSp=repmat([0:L]',(L+1),1);
  vzl=repmat(NaN,1,length(TH));
  
  for inx=1:length(TH)
    for e=0:2*L
      eB(e+1,1)=sum([2*bigS'+1].*[2*bigSp'+1].*...
		    zeroj(bigS,e,bigSp,Leff,[],C0,S0).^2);
    end
    % For the single/double cap
    [BeL(inx,:),dels]=bpboxcap(TH(inx),2*L,[],0,sord);
    vzl(inx)=1/4/pi/K(inx)^2*[(2*[0:2*L]+1).*BeL(inx,:)]*eB;
  end
  % Compare just for good measure
  difer(vzl'-v(:,1))
  disp('Check for zero-l passed')
  % Another check for the zero-l
  for inx=1:length(TH)
    [G,V,EL,EM,KN]=glmalpha(90*(sord==2)+(-1)^(sord+1)*TH(inx),L,sord);
    Vzchk(inx)=sum(V.^2)/sum(V)^2;
  end
  difer(Vzchk'-v(:,1))
  disp('Another check for zero-l passed')
end

% Now do the numerical checks
% What should the large-l limit be? Use the Legendre functions
vll=[1/(4*pi)*[repmat(2*p+1,length(TH),1).*Gp]*plm(0:2:2*L,0,0).^2]';  

yt=sprintf('multitaper variance ratio');
xt=sprintf('degree l');

axes(ah)
% Plot large-l limit
px=plot([0 l(end)],...
	repmat(round(sort(vll)*10^(sord+1))/10^(sord+1),2,1),...
	'Col',grey);
hold on
% Plot large-area limit
%pl=plot(l,vWS,'Col',grey);

% This is where you can explicitly test that the approximation works
Ato4pi=K/(L+1)^2;
% This is the empirical exponent quoted in the paper
empirex=0.88;
plxtra=plot(l,repmat(vWS,length(TH),1)./...
	    repmat(Ato4pi.^empirex,1,length(l)),'Col',grey);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pd=plot(l(end),1./(2*L+1),'k+');
pl=plot(l(end),vWS(end),'o','MarkerS',3,'MarkerE','k','MarkerF','w');

% Plot actual variance ratio
p=plot(l,v,'k','Linew',1);

xmax=l(end);
ymin=0;
% New ymax
ymax=round(max(v(:)*1.15)*10^(sord+1))/10^(sord+1);
ytix=sort([1 ymin ymax]);
ytix=unique([round(vll*10^(sord+1))/10^(sord+1) ymax ymin]);
ylims=[ymin ymax];
ylabs=num2cell(ytix);
% Complete faking
if ylabs{3}==0.5
  ylabs{3}='0.50';
elseif ylabs{end}==0.22
  ylabs{end}='0.220';
end
%ylabs{1}=sprintf('1/(L+1)%s',str2mat(178));

pp=plot([L L],[0 ymax],'k:');

set(ah,'Ytick',ytix,'Ylim',ylims,'ytickl',ylabs,'xlim',[0 xmax])
longticks(ah(1))

set(p,'marker','o','markers',2,'lines','none',...
	     'markerfacec','k','markere','k','lines','-','linew',0.5)

yl=ylabel(yt);
xl=xlabel(xt);

[b,t]=boxtex('ur',ah(1),sprintf('L = %i',L),12,1,0.85,1.15);

% Produce right labels etc - need symbol font, keep it for last
[THsort,i]=sort(TH,'descend'); K=K(i);
for in=1:length(THsort)
  funlob{in}=sprintf('%s = %i%s  K = %i','\Theta',...
		     90*(sord==2)+(-1)^(sord+1)*THsort(in),...
		     str2mat(176),round(K(in)));
  %  tx(in)=text(34.5-(sord==2),v(i(in),end)+(ymax-ymin)/20,funlob{in});
  tx(in)=text(33.5,v(i(in),end)+(ymax-ymin)/20,funlob{in});
end
set(tx,'HorizontalA','left','FontS',8)

set(ah,'box','on')
lah=caplogo(ah,sord+(sord==2),'ur');
movev(lah,-0.053)
moveh(lah,-0.075)

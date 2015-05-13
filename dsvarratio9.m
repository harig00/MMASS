function dsvarratio9
% DSVARRATIO9
%
% Plots the large-l MT variance ratio on a log scale as a function of the
% angle and the degree

clf
[ah,ha]=krijetem(subnum(2,2));

% First two panels show variance ratio in function of degree
% First panel
axes(ah(1))
% Calculate
sord=1;
EL=[0:20];
% This is the cap
TH=[0 10  20 30 40 50 70 100 180];
[SvllMT,K,SvllMTK1,KK1]=doit(EL,TH,sord);
yls=[-1.7 1];
[pl,pd,pp,xl(1),yl(1),tp,pbf1]=...
    plotit(SvllMT,K,SvllMTK1,KK1,EL,TH,sord,yls);
set(pp(1),'LineS','-')
set(pp(2),'LineS','-')

% Second panel
axes(ah(2))
% Calculate
sord=2;
EL=[0:20];
% This is the cut with the equivalent cap size in brackets
TH=90-[0 5 10:10:40 60 90];
[SvllMT,K,SvllMTK2,KK2]=doit(EL,TH,sord);
yls=[-1.7 1];
[pl,pd,pp,xl(2),yl(2),tp,pbf2]=...
   plotit(SvllMT,K,SvllMTK2,KK2,EL,TH,sord,yls);
set(pp(1),'LineS','-')
set(pp(2),'LineS','-')

% Last two panels show variance ratio in function of degree
% Third panel
axes(ah(3))
% Calculate
sord=1;
EL=[0:1:3 5 7 10 14 20];
% This is the cap
TH=1:180;
% Make this go 
[SvllMT,K]=doit(EL,TH,sord);
yls=[-1.7 0.3];
[pl,pd,pp,xl(3),yl(3),tp,ax(1),ps1,pf1,pbfbf1]=...
    plotthat(SvllMT,K,SvllMTK1,KK1,EL,TH,sord,yls);
set(pp(1),'LineS','-')
set(pp(2),'LineS','-')
% Put on an extra axis tick marks at 1/(2L+1)
set(ah(3),'ytick',sort(log10(1./(2*EL+1))),'ygrid','on')
% And see to it that indeed, this is about the value reached by A/4pi
% equals one

% Fourth panel
axes(ah(4))
% Calculate
sord=2;
EL=[0:1:3 5 7 10 14 20];
% This is the cut with the equivalent cap size in brackets
TH=90-[1:90];
[SvllMT,K]=doit(EL,TH,sord);
yls=[-1.7 0.3];
[pl,pd,pp,xl(4),yl(4),tp,ax(2),ps2,pf2,pbfbf2]=...
    plotthat(SvllMT,K,SvllMTK2,KK2,EL,TH,sord,yls);
set(pp(1),'LineS','-')
set(pp(2),'LineS','-')
set(ah(4),'ytick',sort(log10(1./(2*EL+1))),'ygrid','on')

% Cosmetics
delete(yl([2 4]))
nolabels(ah([2 4]),2)
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ax(1:2),1/2,'across')

movev(xl(3:4),-.1)
longticks(ah)

% Put caplogos on
lah(1)=caplogo(ah(1),1,'ll');
lah(2)=caplogo(ah(2),3,'ll');
lah(3)=caplogo(ah(3),1,'ll');
lah(4)=caplogo(ah(4),3,'ll');

set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')

% What's shown here uses the 0.88 rule
delete([ps1(~~ps1) ps2(~~ps2)])
delete([pbfbf1 pbfbf2])
delete([pbf1(:) ; pbf2(:)])

% See also DSVARRATIO

% Show the fits
[pf1 log10(2*EL'+1)];
[pf2 log10(2*EL'+1)];

% Get rid of the superfluous stuff
delete(ah(3:4))
delete(ax)
delete(lah(3:4))

figdisp([],[],[],0)
!degs /home/fjsimons/EPS/dsvarratio9.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SvllMT,K,SvllMTK,KK]=doit(EL,TH,sord)
% This here performs the necessary calculations

% Initialize variance
SvllMT=nan(length(EL),length(TH));

% Load the database just once
[jk,C0,S0]=zeroj(0,0,0,2*max(EL));

% Run this calculation
K=nan(length(EL),length(TH));

% Sort-of axis labels in terms of Shannon numbers
KK=[1 10 100];

indx=0;
for L=EL
  indx=indx+1;
  
  % Calculate the others knowing you must supply the zero-TH one
  [Gpall,p,K(indx,:)]=gammap(L,TH,sord);
  
  if nargout>2
    % Calculate the meaningful Shannon numbers for the labeling, later
    GpK=gammap(L,KK,sord,1,1,1);
  end
  
  % Need zero-area result explicitly provided to you
  if sord==1
    needzero=TH(1)==0 & all(isnan(Gpall(1,:)));
  elseif sord==2
    needzero=TH(end)==0 & all(isnan(Gpall(1,:)));
  end
  if needzero
    % Must put in the TH=0 separately
    % Make the fixed-L null-sphere (TH0) approximation
    bigS=gamini([0:L],(L+1))';
    bigSp=repmat([0:L]',(L+1),1);
    GpTH0=repmat(NaN,1,2*L);
    % This is the curly gamma of yore
    for p=0:2:2*L
      % Watch the square
      GpTH0(p+1)=[sum((2*bigS'+1).*(2*bigSp'+1).*...
		      zeroj(bigS,p,bigSp,2*max(EL),[],C0,S0).^2)]^2;
    end
    % This must be exact 
    GpTH0=GpTH0(1:2:end)*4*pi/(L+1)^4;
    
    % Check that at p=0 this is alwas 4*pi
    difer(GpTH0(1)-4*pi)
    if sord==1
      % Fit the others in there
      Gpall(1,:)=GpTH0;
    elseif sord==2
      % Fit the others in there
      Gpall(1,:)=GpTH0;
    end
  end

  % The even degrees of GAMMAP that matter for the variance
  p=0:2:2*L;
  % Large-l limit of the MT variance ratio
  SvllMT(indx,:)=1/(4*pi)*...
      repmat(2*p+1,length(TH),1).*Gpall*plm(p,0,0).^2;
  if nargout>2
    % Same but for the constant Shannon number
    SvllMTK(indx,:)=1/(4*pi)*...
	repmat(2*p+1,length(KK),1).*GpK*plm(p,0,0).^2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,pd,pp,xl,yl,tp,pbfbf]=...
    plotit(SvllMT,K,SvllMTK,KK,EL,TH,sord,yls)

xel=log10(2*EL+1);
xelm=log10(2*max(EL)+1);

% Plot the variance in function of the Shannon number
% But not the zero-L result since it's ugly
pp=plot(xel(2:end),log10(SvllMTK(2:end,:)),...
	'-','Color',grey,'LineW',1);
hold on

% Plot the actual results
pl=plot(xel,log10(SvllMT),'k','LineW',0.5);

% The last line, for the whole sphere, I believe, should be
%% pp1=plot(xel,log10(1./(2*EL+1)),'r')
% And it's slightly surprising that it isn't
%% pp2=plot(xel,log10(1./K(:,end)),'r')
% I suppose the comparison should be to the whole-sphere estimate after
% that's been binned into bins of (2*L+1) width - after that, having more
% than one taper should reduce that variance by 1/K if they are all
% approximately uncorrelated.

% Put the TH text on
for indx=1:length(TH)
  funlob{indx}=sprintf('%s = %i%s','\Theta',...
		     90*(sord==2)+(-1)^(sord+1)*TH(indx),...
		     str2mat(176));
  tp(indx)=text(xelm+xelm/30,...
		log10(SvllMT(end,indx)),funlob{indx});
end
set(tp,'HorizontalA','left','FontS',8)
% Put the K text on
%if sord==1
  % Effectively delete the first point
  SvllMTK(1,:)=NaN;
  offs=[0 0.1 0.05];
  offs=[-0.2 -0.2 -0.2];
  offx=[-0.2 -0.1 -0.02];
  for indx=1:length(KK)
    funlobk{indx}=sprintf('%s = %i%s','K',KK(indx));
%    tpk(indx)=text(xelm+xelm/30,...
%	 	log10(SvllMTK(end,indx))+offs(indx),funlobk{indx});
    tpk(indx)=text(offx(indx)+xel(1+sum(isnan(SvllMTK(:,indx)))),...
	 	log10(SvllMTK(1+sum(isnan(SvllMTK(:,indx))),indx))...
		   +offs(indx),funlobk{indx});
  end
  set(tpk,'HorizontalA','right','FontS',8,'Color','k')
%end

% Hold on and give some sort of color coding for the Shannon number
ELS=[repmat(EL,length(TH),1)]';
intv=[0 1 10 100 3*max(K(:))];
symbs={'o' 'o' 's' 's'};
for ind=1:length(intv)-1
  selx=K>=intv(ind) & K<intv(ind+1);
  if sum(~~selx(:))
    pt(ind)=plot(log10(2*ELS(selx)+1),log10(SvllMT(selx)),...
		 symbs{ind},'Color','k');
  end
end

if length(pt)>=4
  % Color patches
  set(pt([1 3]),'MarkerF','w','MarkerE','k')
  set(pt([2 4]),'MarkerF','k','MarkerE','k')
  set(pt,'MarkerS',3)
end

% Plot the origin as a special point
pd=plot(0,0,'d','markers',3,'markerf','w','markere','k');

% Plot the best-fit best fit with the area etc
Ao4p=K./repmat((EL'+1).^2,1,size(K,2));
bfbfb=-repmat(xel',1,size(K,2))-0.88*log10(Ao4p);
bfbfb(K<20)=NaN;
pbfbf=plot(xel,bfbfb,'b+');

hold off

xl=xlabel('bandwidth L');
yl=ylabel('variance ratio (\sigma_\infty^2)^{MT}');

% X-axis labeling
% Should write the following into a little function...
% fill cell with string from nums except some
for indx=1:length(EL)
  elab{indx}=num2str(EL(indx));
end
nix=[7:10 12:15 16 17 18 19 20];
for indx=1:length(nix)
  elab{nix(indx)}='';
end
ontop=0;
if ontop
  [ax,xl,yl]=xtraxis(gca,xel,elab,'L');
  longticks(ax)
else
  set(gca,'xtick',xel,'xtickl',elab)
end

% Y-axis labeling
syls1=[0.01:0.01:0.1];
syls2=[0.1:0.1:1];
syls3=[1:1:10];
syls=unique([syls1 syls2 syls3]);
for indx=1:length(syls)
  elab{indx}=num2str(syls(indx));
end
nix=[2:9 11:18 20:27];
for indx=1:length(nix)
  elab{nix(indx)}='';
end
set(gca,'ytick',log10(syls),'ytickl',elab)

ylim(yls)
xlim([0 xelm+xelm/3.5])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,pd,pp,xl,yl,tp,ax,ps,pf,pbfbf]=...
    plotthat(SvllMT,K,SvllMTK,KK,EL,TH,sord,yls)

% Figure out the area
A=4*pi*spharea(TH,sord);
if sord==3
  error('This choice is not allowed')
end

xel=log10(A/4/pi);
xelm=0;
xelmi=-2;
AK=KK(:)*[4*pi./([1:size(SvllMTK,1)-1]+1).^2];
% Plot the variance in function of the Shannon number
pp=plot(log10(AK/4/pi)',log10(SvllMTK(2:end,:)),'-',...
	'Color',grey,'LineW',1);
hold on

% Plot the actual results
pl=plot(xel,log10(SvllMT),'k','LineW',0.5);

% Fit the lines through now whenever the K's are bigger than 10
for indx=1:size(SvllMT,1)
  Ksel=K(indx,:)>=10;
  if ~~sum(Ksel)
    % Row index into SvllMT is for different bandwidths
    pf(indx,:)=polyfit(xel(Ksel),log10(SvllMT(indx,Ksel)),1);
    hold on
    % Plot the prediction... pf(:,1) is the slope, pf(:,2) the intercept
    ps(indx)=plot(xel(Ksel),...
		  pf(indx,2)+pf(indx,1)*xel(Ksel),'b-','linew',2);
  end
end

% Plot the best-fit best fit
bfbfb=-log10(repmat(2*EL'+1,1,length(xel)))...
	   -0.88*repmat(xel,length(EL),1);
bfbfb(K<20)=NaN;
pbfbf=plot(xel,bfbfb,'b+');

% Put the text on
for indx=1:length(EL)
  funlob{indx}=sprintf('L = %i',EL(indx));
  if sord==1
    wats=log10(SvllMT(indx,end));
  elseif sord==2
    wats=log10(SvllMT(indx,end));
  end
    tp(indx)=text(xelm+abs(xelmi)/30,wats,funlob{indx});
end
set(tp,'HorizontalA','left','FontS',8)

plotdots=0;
if plotdots
  % % Hold on and give some sort of color coding for the Shannon number
  AS=repmat(A/4/pi,length(EL),1);
  intv=[0 1 10 100 3*max(K(:))];
  symbs={'o' 'o' 's' 's'};
  for ind=1:length(intv)-1
    selx=K>=intv(ind) & K<intv(ind+1);
    if sum(~~selx(:))
      pt(ind)=plot(log10(AS(selx)),log10(SvllMT(selx)),...
  		 symbs{ind},'Color','k');
    end
  end
  if length(pt)>=4
    % Color patches
    set(pt([1 3]),'MarkerF','w','MarkerE','k')
    set(pt([2 4]),'MarkerF','k','MarkerE','k')
    set(pt,'MarkerS',3)
  end
end

hold off

xl=xlabel('fractional area A/(4\pi)');
yl=ylabel('variance ratio (\sigma_\infty^2)^{MT}');

% X-axis labeling
sxls1=[0.001:0.001:0.01];
sxls2=[0.01:0.01:0.1];
sxls3=[0.1:0.1:1];
sxls=unique([sxls1 sxls2 sxls3]);
for indx=1:length(sxls)
  elab{indx}=num2str(sxls(indx));
end
nix=[2:9 11:18 20:27];
for indx=1:length(nix)
  elab{nix(indx)}='';
end
set(gca,'xtick',log10(sxls),'xtickl',elab)

% Y-axis labeling
syls1=[0.01:0.01:0.1];
syls2=[0.1:0.1:1];
syls3=[1:1:10];
syls=unique([syls1 syls2 syls3]);
for indx=1:length(syls)
  elab{indx}=num2str(syls(indx));
end
nix=[2:9 11:18 20:27];
for indx=1:length(nix)
  elab{nix(indx)}='';
end
set(gca,'ytick',log10(syls),'ytickl',elab)

xlim([xelmi xelm+abs(xelmi)/4.5])
ylim(yls)

switch sord
 case 1
  % This is the cap size
  xTH=[1 10:10:50 70 90 180];
  A=2*pi*(1-cos(xTH*pi/180));
 case 2
  % This is the cut size with the cap size in brackets
  xTH=90-[1 10:10:40 60 90];
  A=4*pi*(1-sin(xTH*pi/180));
 otherwise
  error('This choice is not allowed')
end
% Here want to put a TH axis as well
for indx=1:length(xTH)
  elab{indx}=sprintf('%i%s',...
		     (90*(sord==2)+(-1)^(sord+1)*xTH(indx)),...
		     str2mat(176));
end
[ax,xl2,yl2]=xtraxis(gca,log10(A/4/pi),elab,...
		   sprintf('%s','\Theta'));
longticks(ax)
set(ax,'FontS',8)
movev(xl2,-range(ylim)/7.75)
moveh(xl2,range(xlim)/2.2)

pd=NaN;

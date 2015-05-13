function dsvarratio(l)
% DSVARRATIO(l)
%
% Plots the variance ratio of the periodogram to the whole-sky estimate
%
% INPUT:
%
% l       degrees at which the variance ratio is sought (vector)
%
% Last modified by fjsimons-at-alum.mit.edu, 4/29/2007

% Lmax is automatically fixated in PERIODOVAR

defval('l',0:50)
defval('meth',1)
defval('Lmax',2*max(l))
% The single/double cap are exactly the same - prove it or don't
% 'cos the area is double, but the coupling is over half the number of
% degrees 
defval('xver',0)

% Size of the spherical caps
TH=[3 4 5 7 10 20 60];

% Initialize arrays
v=repmat(NaN,length(l),length(TH));

% SINGLE CAP
% First time return the spectra as well - if you specified Lmax
 [v(1,:),vll,vsth,vlth,TH,A,Bl,bls]=periodovar(l(1),TH,[],[],1,[],Lmax);
for index=2:length(l)
  v(index,:)=periodovar(l(index),TH,Bl,bls,1,[],Lmax);
end
% % Find appropriate Lmax? Hope you've reached the point where Bl is tiny
% for in=1:size(Bl,2)
%   sk=Bl(:,in)<max(Bl(:,in))/1000;
%   if sum(sk)
%     LmaxN(in)=indeks(find(sk),1);
%   else
%     LmaxN(in)=NaN;
%   end
% end
% if any(Lmax<LmaxN)
%   error('Reconsider and make Lmax bigger')
% end

% Do the plotting
clf

ytag1=sprintf('periodogram variance ratio');
ytagx=sprintf('single cap radius %s','\Theta');
ytagx2=sprintf('double cap radius %s','\Theta');
ytags=sprintf('single or double cap radius %s','\Theta');
xtag=sprintf('degree l');
ylima=[0 25];

[ah,ha]=krijetem(subnum(2,2));

axes(ah(1))
px=plot([0 l(end)],repmat(round(sort(vll)*10)/10,2,1),'Col',grey);
hold on
p(1)=plot(l,2*l+1,'LineS','-','Color',grey);
p=[p ; plot(l,v,'k','Linew',1)];

ms=1.5;

plot(0,1,'o','MarkerF','w','MarkerE','k','markers',ms*2)

ylim(ylima)
delete(ah(3:4)); ah=ah(1:2);
movev(ah,-0.2)
yl(1)=ylabel(ytag1);
xl(1)=xlabel(xtag);
set(ah(1),'box','off','ytick',[1 5:5:ylima(2)])
funlab=sprintf('%s = %i%s',' Q',0,str2mat(176));
[THsort,i]=sort(TH,'descend');
for in=1:length(THsort)
  funlob{in}=sprintf('%s = %i%s',' Q',THsort(in),str2mat(176));
end
[ax,xlx,ylx]=xtraxis(ah(1),(ylima(2)-1)/2,funlab,...
		     '',sort(v(end,:)),funlob,ytagx);
set(ax,'fontn','symbol')
set(ylx,'rotation',-90)
% Set Y tick marks to the large-l limit
set(ah(1),'ytick',[round(sort(vll)*10)/10 ylima(2)],'ygrid','off')

if xver==1
  lah=caplogo(ah(1),1,'ul');
  movev(lah,-0.055)
  % Now for the double cap ... this is exactly the same!
  axes(ah(2))
  % DOUBLE CAP
  % First time return the spectra as well... it's about how much is subtracted
  [v(1,:),vll,vsth,vlth,TH,A,Bl,bls]=periodovar(l(1),90-TH,[],[],2,[],Lmax);
  % Watch out! The TH is being returned now
  for index=2:length(l)
    v(index,:)=periodovar(l(index),TH,Bl,bls,2,[],Lmax);
  end
%   % Find appropriate Lmax? Hope you've reached the point where Bl is tiny
%   for in=1:size(Bl,2)
%     sk=Bl(:,in)<max(Bl(:,in))/1000;
%     if sum(sk)
%       LmaxN(in)=indeks(find(sk),1);
%     else
%       LmaxN(in)=NaN;
%     end
%   end
%   if any(Lmax<LmaxN)
%     error('Reconsider and make Lmax bigger')
%   end
  p(end+1)=plot(l,2*l+1,'LineS','-','Color',grey);
  hold on
  p=plot(l,v,'k');
  set(p,'LineW',1)
  plot(0,1,'o','MarkerF','w','MarkerE','k','markers',ms)
  
  ylim(ylima)
  yl(2)=ylabel(ytag1);
  xl(2)=xlabel(xtag);
  set(ah(1),'box','off','ytick',[1 5:5:ylima(2)])
  [ax2,xlx2,ylx2]=xtraxis(ah(2),(ylima(2)-1)/2,0,'',sort(v(end,:)),...
		       sort(90-TH,'descend'),ytagx2);
  set(ylx2,'rotation',-90)
  
  % Set Y tick marks to the large-l limit
  set(ah(2),'ytick',[round(sort(vll)*10)/10 ylima(2)],'ygrid','off')
  lah=caplogo(ah(2),3,'ul');
  movev(lah,-0.055)
  nolabels(ah(2),2)
  delete(yl(2))
  longticks([ah ax2])
  moveh([ah(1) ax],.01)
  deggies(ax2,2)
else
  delete(ah(2))
  ah=ah(1);
  moveh([ah(1) ax] ,.2)
  delete(ylx)
%  set(ylx,'string',ytags)
end

% This needs to come on the bottom
set(p(2:end),'marker','o','markers',ms,'lines','-',...
	     'linew',0.5,'markerfacec','k','markere','k')

% Cosmetics applicable to both
longticks(ah)
set(ax,'TickDir','out','TickLength',[0 0])

% Get rid of double axis and provide alternative
delete(ax)

% Calculate the fraction of the sky that survives the cut
Acap=2*pi*(1-cos(TH*pi/180));

% Do our own minimization here
G=[repmat(1,length(Acap),1) sqrt(4*pi./Acap)'];
svown=inv(G'*G)*G'*vll';
% And compare with polyfit
% sv=polyfit(sqrt(4*pi./Acap),vll,1)
% The first is the slope ; sval=sv(1);
% But really we want a constrained fit through (0,0)
F=[1 0];
d=[G'*vll' ; 0];
G=[G'*G F' ; F 0];
svown=inv(G)*d;
% The first two are the model parameters, the third is the useless
% Lagrange parameter; the second is the slope, the first the zero
% intercept
sval=svown(2);

% Newest thing: add sky fraction
[ax,xlx,ylx]=xtraxis(ah(1),[],[],[],...
		     sval.*sqrt(4*pi./Acap),...
		     round(sval.*sqrt(4*pi./Acap)*10)/10,...
		     sprintf('%3.2f (4%s/A^{cap})^{1/2}',sval,'\pi'));
set(ylx,'rotation',-90)
moveh(ylx,3)
longticks(ax)

%set(ha(1),'box','on')
sord=1;
for in=1:length(THsort)
  funlob{in}=sprintf('%s = %i%s','\Theta',...
		     90*(sord==2)+(-1)^(sord+1)*THsort(in),...
		     str2mat(176));
  tx(in)=text(43.5,v(end,i(in))+range(ylima)/25,funlob{in});
end
set([xl yl xlx ylx],'FontS',8)
set([ah ax],'FontS',7)
set(tx,'HorizontalA','left','FontS',6)

set(gcf,'color','w','inverthardcopy','off')
fig2print(gcf,'portrait')
figdisp([],[],[],0)
!degs /home/fjsimons/EPS/dsvarratio.eps

disp('Make sure the top label is right - fix, or rerun if necessary')

function varargout=dsvarratio2(Lmax,TH,invm)
% [v,TH]=DSVARRATIO2(Lmax,TH,invm)
%
% Plots the variance ratio of the debiased periodogram to the whole-sky
% estimate - for DOUBLE caps only (if not, cannot invert coupling matrix)
%
% INPUT:
%
% Lmax     the maximal degree of the boxcar coupling matrix
% TH       half the belt size in degrees
% invm     inversion method, a string
%          'inv' [default], 'pinv' or 'geninv' 
%          all these methods agree for low condition numbers
%
% OUTPUT:
%
% v        The variance ratios
% TH       half the belt size used, in degrees
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2007

% This is the inverse of the symmetric part of the boxcar coupling matrix

% Double cap is the only relevant computable quantity
sord=2;
defval('Lmax',150);
defval('TH',[1  5 10 15]);
defval('xver',0)
defval('invm','inv')
ylima=[0.9 1.99];
xlima=[0 50];

% Calculate the fraction of the sky that survives the cut
fsky=(1-sin(TH*pi/180));

% If the condition number exceeds log10(eps)... it's rubbish ...
% but actually, this seems to happen much earlier... keep it near one 
v=repmat(NaN,Lmax+1,length(TH));
for index=1:length(TH)
  [K,A]=bcoupling(TH(index),Lmax,sord);
  disp(sprintf('condition number %4.2f',cond(K)))
  % Now invert this guy, find diagonal, multiply by 4pi and divide by 2l+1 
  v(:,index)=diag(feval(invm,K))*4*pi/A;
  if xver==1
    % Try this other way - same deal
    K=K./repmat(2*[0:Lmax]+1,Lmax+1,1)*A;
    % I think this condition number is always higher
    disp(sprintf('modified condition number %4.2f',cond(K)))
    v2=diag(feval(invm,K))*4*pi./(2*[0:Lmax]'+1);
    difer(v(index,:)-v2)
  end
end

clf
[ah,ha]=krijetem(subnum(2,2));

axes(ah(1))

p=plot(0:Lmax,v,'k');
ylim(ylima)
xlim(xlima)
hold on
plot(xlima,[1 1],'Color',grey)

% Cosmetics
delete(ah(2:4)); ah=ah(1);
moveh(ah,.2)
movev(ah,-0.2)
[THsort,i]=sort(90-TH,'descend');
for in=1:length(THsort)
  funlob{in}=sprintf('%s = %i%s',' Q',THsort(in),str2mat(176));
end
xtag=sprintf('degree l');
ytag1=sprintf('maximum likelihood variance ratio');
yl(1)=ylabel(ytag1);
xl(1)=xlabel(xtag);

[ax,xlx,ylx]=xtraxis(ah(1),[],[],[],v(xlima(2)+1,:),funlob);
set(ax,'fontn','symbol')
set(ylx,'rotation',-90)

ms=1.5;

set(p,'marker','o','markers',ms,'lines','-',...
	     'linew',0.5,'markerfacec','k','markere','k')
longticks(ah)

% New cosmetics
delete(ax)

% Newest thing: add sky fraction
[ax,xlx,ylx]=xtraxis(ah(1),[],[],[],...
		     1./(fsky.^2),round(1./(fsky.^2)*100)/100,...
		     sprintf('(4%s/A)^2','\pi'));
set(ylx,'rotation',-90)
moveh(ylx,3)
longticks(ax)

%set(ha(1),'box','on')
for in=1:length(THsort)
  funlob{in}=sprintf('%s = %i%s','\Theta',...
		     THsort(in),...
		     str2mat(176));
  tx(in)=text(43.5,v(end,i(in))+range(ylima)/20,funlob{in});
end
set([xl yl xlx ylx],'FontS',8)
set([ah ax],'FontS',7)
set(tx,'HorizontalA','left','FontS',6)

% Put on caplogo last
lah=caplogo(ah(1),3,'ul');
movev(lah,-0.055)

set(gcf,'color','w','inverthardcopy','on')
figdisp([],[],[],0)
!degs /home/fjsimons/EPS/dsvarratio2.eps

varns={v,TH};
varargout=varns(1:nargout);

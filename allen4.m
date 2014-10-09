function allen4(sels)
% ALLEN4(sels)
%
% Plots the wavelet diagnostics for early warning
%
% INPUT:
%
% sels        0 No distance selection
%             1 Distance selection to within 100 km
%
% Last modified by fjsimons-at-alum.mit.edu, 11/30/2006

defval('sels',0)

stuff=load('/home/fjsimons/EALARMS/wtALL_diagnostics');

% Collect the data
mom=stuff(:,4);
dist=stuff(:,5);
evnu(:,1)=stuff(:,6);
skal(:,1)=stuff(:,7);
skal(:,2)=stuff(:,8);
skal(:,3)=stuff(:,9);
skal(:,4)=stuff(:,10);
skal(:,5)=stuff(:,11);

% Protect identical moments belonging to individual events
% Don't do this and you get what's in the paper
% Do this and you work per event, not per magnitude
% Do this and you have to adjust the breakpoint
% mom=mom+evnu/1000;

for index=1:5
  nansf(:,index)=isnan(skal(:,index));
  nuna(index)=round(100*sum(nansf(:,index))/length(mom));
  disp(sprintf('Number of NaNs at scale %2.2i is %i%s',...
	       index,nuna(index),'%'))
  % Try to figure out where we are missing the data...
  % Protect against feeding logicals
  [H,c11,cmn]=bindens(mom,0+nansf(:,index),[],1);
  watm(:,index)=H(1,:)*100;
end
% Later can do 
% plot(linspace(c11(1),cmn(1),size(watm,1)),watm)
% The number of totally missed events
disp(sprintf('The number of missed events is %i',...
	     sum(sum(isnan(skal),2)==5)))

% First four or 2 through 5?
offs=1;

% Thresholding division
thresh=stuff(:,3);
difer(diff(thresh))
disp(sprintf('Thresholding divisor %2.2f',thresh(1,1)))

if sels==1
  % Take out those outside of 100 km
  selkt=dist<100;
  mom=mom(selkt);
  skal=skal(selkt,:);
end

clf
[ah,ha]=krijetem(subnum(2,2));
txt1=[5.1 60];
txt2=[5.1 17.5];
txt3=[5.1 220];
txt4=[6.9 1000];

% This is the magnitude breakpoint for the regression
slev=20;

for index=[1:4]+offs
  % Plot the coefficients in skal1 already logarithmically scaled
  % Get the statistics
  [s,g]=row2kstats(mom,log10(skal(:,index)));
  axes(ah(index-offs))
  % Do this for visual only...
  % nans=isnan(skal(:,index));
  % [H,c11,cmn]=bindens(mom(~nans),log10(skal(~nans,index)));
  % imagef(c11,cmn,H)
  % Plot the individual grey dots for the actual data
  t(:,index)=semilogy(mom,skal(:,index),'o');
  hold on
  % Plot mean and median; get rid of rest
  [p(:,index),u(:,index),jk1,jk2]=plotsingle(s,g);
  delete([jk1 jk2])
  % Plot box plot; get rid of rest
  [jk1,jk2,m(:,index),l(:,index)]=plotsingle(s,g);
  delete([jk1 jk2])
  % Now the fit to the data for the low magnitudes
  Gmeanl=polyfit(s(1:slev),g.mean(1:slev),1);
  % Magnitude break point
  disp(sprintf('Magnitude breakpoint %4.2f',s(slev)))
  yfit=polyval(Gmeanl,s(1:slev));
  % Calculate the correlation coefficient
  [Rl,pl]=corrcoef(s(1:slev),g.mean(1:slev));
  % Plot the regression line
  fpl(index)=semilogy(s(1:slev),10.^yfit,'k-','LineW',1);
  % Now the fit to the data for the high magnitudes
  Gmeanh=polyfit(s(slev+1:end),g.mean(slev+1:end),1);
  yfit=polyval(Gmeanh,s(slev+1:end));
  % Calculate the correlation coefficient
  [Rh,ph]=corrcoef(s(slev+1:end),g.mean(slev+1:end));
  % Plot the regression line
  fph(index)=semilogy(s(slev+1:end),10.^yfit,'k-','LineW',1);
  % Annotate plot
  tl(index)=text(txt1(1),txt1(2),...
		 sprintf('m_l = %4.2f logC_%i %+4.1f',...
			 1/Gmeanl(1),index,...
			 -Gmeanl(2)/Gmeanl(1)));
  th(index)=text(txt2(1),txt2(2),...
		 sprintf('m_h= %4.2f logC_%i %+4.1f',...
			 1/Gmeanh(1),index,...
			 -Gmeanh(2)/Gmeanh(1)));
  % 5 % confidence level for significance
  if pl(2)>0.05; Rl(2)=NaN; end
  if ph(2)>0.05; Rh(2)=NaN; end  
  tr(index)=text(txt3(1),txt3(2),...
		 sprintf('R_l = %4.2f    R_h= %4.2f',...
			 Rl(2),Rh(2))); 
  tr(index)=text(txt4(1),txt4(2),...
		 sprintf('%i%s',100-nuna(index),'%')); 
  hold off
  axis tight
  xl(index-offs)=xlabel('TriNet magnitude');
  yl(index-offs)=ylabel('wavelet coefficient magnitude');  
end

% Cosmetics
set(p(~~p),'MarkerF','k','MarkerS',3)
set(u(~~u),'MarkerF','y','MarkerS',3)
set(t(~~t),'MarkerF',grey,'MarkerS',2,'MarkerE',grey)
set(ah,'xlim',[2.8 7.5],'ylim',[10 10^6.75])
set(ah,'ytick',10.^[2:2:6])
nolabels(ah(1:2),1)
nolabels(ha(3:4),2)
delete(xl(1:2))
delete(yl([2 4]))
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')
longticks(ah)
delete(u(~~u))
fig2print(gcf,'portrait')
% [bh,th]=label(ah,'lr',12,[],2);
for index=1:length(ah)
  % LABELS ON LOGARITHMIC PLOTS ARE A PAIN
  % INCOROPRATE THIS INTO THE LABEL ROUTINE
  [ax(index),axl(index,:)]=laxis(ah(index),0);
  xll(index,:)=get(ax(index),'xlim');
  yll(index,:)=get(ax(index),'ylim');
  [bh,th]=label(ax,'ul',12,offs,2);
  set(ax(index),'xlim',xll(index,:));
  set(ax(index),'ylim',yll(index,:));
end

figdisp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pp,uu,mm,ll]=plotsingle(ss,gg,limz)
defval('limz',length(ss))
pp=semilogy(ss(1:limz),10.^gg.mean(1:limz),'ko');
uu=semilogy(ss(1:limz),10.^gg.median(1:limz),'ko');
levl=0.05;
mm=semilogy(gamini(ss(1:limz),4)+repmat(levl*[-1 0 1 NaN],1,limz),...
	    gamini(10.^gg.p25(1:limz),4),'k-');
ll=semilogy(gamini(ss(1:limz),4)+repmat(levl*[-1 0 1 NaN],1,limz),...
	    gamini(10.^gg.p75(1:limz),4),'k-');


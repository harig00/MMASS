function allen5(sels)
% ALLEN5(sels)
%
% Are the errors uncorrelated?
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

% First four or 2 to 5?
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
slev=20;
% Just calculate the magnitude errors
for index=[1:4]+offs
  axes(ah(index-offs))

  % Get the statistics on the logarithmically scaled data
  [s,g]=row2stats(mom,log10(skal(:,index)));
  % Now the fit to the data for the low magnitudes
  Gmeanl=polyfit(s(1:slev),g.mean(1:slev),1);
  % Magnitude break point
  disp(sprintf('Magnitude breakpoint %4.2f',s(slev)))
  yfit=polyval(Gmeanl,s(1:slev));
  % Figure out the error in the determination
  % Actual seismic magnitude is in s(1:slev)
  % Predicted moment is, from the averaged coefficients and the
  % regression line; predict ALL with the LOW line
  predmoml(:,index)=1/Gmeanl(1)*g.mean-Gmeanl(2)/Gmeanl(1);
  % Plot the error
  % Predicting the LOW magnitudes with the LOW regression
  errll(:,index)=s(1:slev)-predmoml(1:slev,index)';
  % Predicting the HIGH magnitudes with the LOW regression
  errlh(:,index)=s(slev+1:end)-predmoml(slev+1:end,index)';
  ell(index)=plot(s(1:slev),errll(:,index),'s');
  hold on
  elh(index)=plot(s(slev+1:end),errlh(:,index),'s');
  % Now the fit to the data for the high magnitudes
  Gmeanh=polyfit(s(slev+1:end),g.mean(slev+1:end),1);
  yfit=polyval(Gmeanh,s(slev+1:end));
  % Figure out the error in the determination
  % Predict ALL with the high line
  predmomh(:,index)=1/Gmeanh(1)*g.mean-Gmeanh(2)/Gmeanh(1);
  % Plot the error; 
  % Predicting the LOW magnitudes with the HIGH regression
  errhl(:,index)=s(1:slev)-predmomh(1:slev,index)';
  % Predicting the HIGH magnitudes with the HIGH regression
  errhh(:,index)=s(slev+1:end)-predmomh(slev+1:end,index)';
  ehl(index)=plot(s(1:slev),errhl(:,index),'^');
  ehh(index)=plot(s(slev+1:end),errhh(:,index),'^');
  % Labels, etc
  xl(index-offs)=xlabel('TriNet magnitude');
  yl(index-offs)=ylabel('magnitude prediction error');
  grid on
end

% Now figure out the correlation structure between the errors
confid=0.05;
[R,P]=corrcoef([errll(:,1+offs:end) ; ... 
		errlh(:,1+offs:end) ; ... 
		errhl(:,1+offs:end) ; ... 
		errhh(:,1+offs:end)]);
R(P>confid)=NaN;

% This shows that scales 2 and 5 are completely decorrelated

% Cosmetics
set(ehl(~~ehl),'MarkerF',grey,'MarkerS',3,'MarkerE',grey)
set(elh(~~elh),'MarkerF',grey,'MarkerS',3,'MarkerE',grey)
set(ehh(~~ehh),'MarkerF','k','MarkerS',3,'MarkerE','k')
set(ell(~~ell),'MarkerF','k','MarkerS',3,'MarkerE','k')
nolabels(ah(1:2),1)
nolabels(ha(3:4),2)
delete(xl(1:2))
delete(yl([2 4]))
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

set(ah,'ylim',[-4.5 4.5],'ytick',[-4:2:4],'xlim',[2.8 7.5])

longticks(ah)
fig2print(gcf,'portrait')
for index=1:length(ah)
  % LABELS ON LOGARITHMIC PLOTS ARE A PAIN
  % INCOROPRATE THIS INTO THE LABEL ROUTINE
  [ax(index),axl(index,:)]=laxis(ah(index),0);
  xll(index,:)=get(ax(index),'xlim');
  yll(index,:)=get(ax(index),'ylim');
  [bh,th]=label(ax,'ur',12,offs,2);
  set(ax(index),'xlim',xll(index,:));
  set(ax(index),'ylim',yll(index,:));
end

figdisp([])

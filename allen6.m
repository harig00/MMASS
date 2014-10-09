function allen6(sels)
% ALLEN6(sels)
%
% Magnitude prediction from wavelet analysis
% FROM SCALES 2 and 5 AVERAGED, worse than only 5
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
%mom=mom+evnu/1000;

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

slev=20;
% Just calculate the magnitude errors
for index=[1:4]+offs
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
  % Predicting the LOW magnitudes with the LOW regression
  errll(:,index)=s(1:slev)-predmoml(1:slev,index)';
  % Predicting the HIGH magnitudes with the LOW regression
  errlh(:,index)=s(slev+1:end)-predmoml(slev+1:end,index)';
  
  % Now the fit to the data for the high magnitudes
  Gmeanh=polyfit(s(slev+1:end),g.mean(slev+1:end),1);
  yfit=polyval(Gmeanh,s(slev+1:end));
  % Figure out the error in the determination
  % Predict ALL with the high line
  predmomh(:,index)=1/Gmeanh(1)*g.mean-Gmeanh(2)/Gmeanh(1);
  % Predicting the LOW magnitudes with the HIGH regression
  errhl(:,index)=s(1:slev)-predmomh(1:slev,index)';
  % Predicting the HIGH magnitudes with the HIGH regression
  errhh(:,index)=s(slev+1:end)-predmomh(slev+1:end,index)';
end


% Magnitude determination from scales 2 and 5 separately
av25l=mean(predmoml(:,[2 5]),2)';
av25h=mean(predmomh(:,[2 5]),2)';
% Predict LOW magnitudes with LOW average regression
av25errll=s(1:slev)-av25l(1:slev);
% Predict HIGH magnitudes with LOW average regression
av25errlh=s(slev+1:end)-av25l(slev+1:end);
% Predict LOW magnitudes with HIGH average regression
av25errhl=s(1:slev)-av25h(1:slev);
% Predict HIGH magnitudes with LOW average regression
av25errhh=s(slev+1:end)-av25h(slev+1:end);
% Magnitude determination from the AVERAGE of the HIGH and LOW regression
% lines over scales 2 and 5
av25avlh=mean([av25l(:) av25h(:)],2);

plot(s(1:slev),av25errll,'s','Markerf','k','MarkerS',3,'MarkerE','k')

clf
[ah,ha]=krijetem(subnum(3,1));
axes(ah(1))
plot(s(slev+1:end),av25errlh,'s','Markerf',grey,'MarkerS',3,'MarkerE',grey)
hold on
plot(s(1:slev),av25errhl,'^','Markerf',grey,'MarkerS',3,'MarkerE',grey)
plot(s(slev+1:end),av25errhh,'^','Markerf','k','MarkerS',3,'MarkerE','k')

plot(s,s-av25avlh','o','Markerf','k','MarkerS',3,'MarkerE','k')
delete(ah(2:3))

xl(1)=xlabel('TriNet magnitude');
yl(1)=ylabel('magnitude prediction error');
longticks(gca,2)
grid on

set(ah(1),'ylim',[-2 2],'ytick',[-2:2],'xlim',[2.8 7.4])

fig2print(gcf,'portrait')
figdisp



% If it's positive, it's UNDEPREDICTED... PREDICTION IS TOO LOW

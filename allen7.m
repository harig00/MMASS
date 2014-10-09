function allen7(sels)
% ALLEN7(sels)
%
% Magnitude prediction from scale 5 ONLY
% as a function of the TriNet magnitude
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
% Best regression line is when you DO NOT do this
mom=mom+evnu/1e10;

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
disp(sprintf('Work with %i records',length(mom)))

clf
ah=axes;
slev=20;
% Just calculate the magnitude errors at SCALE 5
index=5;

axes(ah)

% Get the statistics on the logarithmically scaled data
[s,g]=row2stats(mom,log10(skal(:,index)));
% Now the fit to the data for the low magnitudes
Gmeanl=polyfit(s(1:slev),g.mean(1:slev),1)
% Use the average regression line on the individual magnitudes instead
Gmeanl=[0.9638   -0.4546]
% Magnitude break point
disp(sprintf('Magnitude breakpoint %4.2f',s(slev)))
yfit=polyval(Gmeanl,s(1:slev));
% Figure out the error in the determination
% Actual seismic magnitude is in s(1:slev)
% Predicted moment is, from the averaged coefficients and the
% regression line; predict ALL with the LOW line
predmoml=1/Gmeanl(1)*g.mean-Gmeanl(2)/Gmeanl(1);
% Make sure this corresponds with what's plotted...
disp(sprintf('m_l = %4.2f logC_%i %+4.1f',...
			 1/Gmeanl(1),index,...
			 -Gmeanl(2)/Gmeanl(1)))
% Plot the error versus the TriNet magnitude
% Predicting the LOW magnitudes with the LOW regression
errll=s(1:slev)-predmoml(1:slev);
% Predicting the HIGH magnitudes with the LOW regression
errlh=s(slev+1:end)-predmoml(slev+1:end);
ell=plot(s(1:slev),errll,'s');
hold on
elh=plot(s(slev+1:end),errlh,'s');
% Now the fit to the data for the high magnitudes
Gmeanh=polyfit(s(slev+1:end),g.mean(slev+1:end),1)
Gmeanh=[0.6851    0.8416]
disp(sprintf('m_h= %4.2f logC_%i %+4.1f',...
			 1/Gmeanh(1),index,...
			 -Gmeanh(2)/Gmeanh(1)))
yfit=polyval(Gmeanh,s(slev+1:end));
% Figure out the error in the determination
% Predict ALL with the high line
predmomh=1/Gmeanh(1)*g.mean-Gmeanh(2)/Gmeanh(1);
% Plot the error; 
% Predicting the LOW magnitudes with the HIGH regression
errhl=s(1:slev)-predmomh(1:slev);
% Predicting the HIGH magnitudes with the HIGH regression
errhh=s(slev+1:end)-predmomh(slev+1:end);
ehl=plot(s(1:slev),errhl,'^');
ehh=plot(s(slev+1:end),errhh,'^');
% Labels, etc
xl=xlabel('TriNet magnitude');
yl=ylabel('magnitude prediction error');
grid on

% Now the average of both regression lines
avhl=mean([predmoml(:) predmomh(:)],2)';
pav=plot(s,s-avhl,'o','Markerf','k','MarkerS',3,'MarkerE','k');

% Cosmetics
hl=2;
set(ehl,'MarkerF',grey,'MarkerS',hl,'MarkerE',grey)
set(elh,'MarkerF',grey,'MarkerS',hl,'MarkerE',grey)
set(ehh,'MarkerF','k','MarkerS',hl,'MarkerE','k')
set(ell,'MarkerF','k','MarkerS',hl,'MarkerE','k')
%delete([ehl elh ehh ell])
%delete(pav)

% This error IS symmetric, it's about their average
[perb,k,l]=errorxy(s,s-avhl,[],abs([errhl errhh]-(s-avhl)));
delete(perb)

shrink(gca,1,3)
set(ah(1),'ylim',[-1.3 1.3],'ytick',[-1:1],'xlim',[2.8 7.4])
longticks(ah,2)
hold off
%g.nonnans

figdisp

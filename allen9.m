function allen9(sels)
% ALLEN9(sels)
%
% Magnitude prediction from scale 5 ONLY
% as a function of the TriNet magnitude,
% the number of stations, and the number of stations versus
% TriNet magnitude
%
% INPUT:
%
% sels        0 No distance selection
%             1 Distance selection to within 100 km
%
% Last modified by fjsimons-at-alum.mit.edu, 18.07.2006

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
% so use the not-done regression but plot them for yes-done events
% We thus define the best regression possible, but show the results
% on individual events, not averaging events with the same magnitude 
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

clf
ah=krijetem(subnum(3,1));
slev=20;
% Just calculate the magnitude errors at SCALE 5
index=5;

% Get the statistics on the logarithmically scaled data
[s,g]=row2stats(mom,log10(skal(:,index)));
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
% Now the fit to the data for the high magnitudes
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
% Now the average of both regression lines
avhl=mean([predmoml(:) predmomh(:)],2)';

% ERROR VERSUS MAGNITUDE
axes(ah(1))
ell=plot(s(1:slev),errll,'s');
hold on
[perb,k,l]=errorxy(s,s-avhl,[],abs([errhl errhh]-(s-avhl)));
hold on
elh=plot(s(slev+1:end),errlh,'s');
ehl=plot(s(1:slev),errhl,'^');
ehh=plot(s(slev+1:end),errhh,'^');
% Labels, etc
xl=xlabel('TriNet magnitude');
yl=ylabel('mag prediction error');
grid on
pav=plot(s,s-avhl,'o','Markerf','k','MarkerS',3,'MarkerE','k');
% This error IS symmetric, it's about their average
hold off
% Cosmetics
hl=2;
set(ehl,'MarkerF',grey,'MarkerS',hl,'MarkerE',grey)
set(elh,'MarkerF',grey,'MarkerS',hl,'MarkerE',grey)
set(ehh,'MarkerF','k','MarkerS',hl,'MarkerE','k')
set(ell,'MarkerF','k','MarkerS',hl,'MarkerE','k')
%delete([ehl elh ehh ell])
%delete(pav)
delete(perb)
set(ah(1),'ylim',[-1.3 1.3],'ytick',[-1:1],'xlim',[2.8 7.4])

% NUMBER OF STATIONS VERSUS TRINET MAGNITUDE
axes(ah(2))
nsvm=plot(s,g.nonnans,'d','MarkerS',hl+1,'MarkerE','k','MarkerF','k');
set(ah(2),'ylim',[-5 70],'ytick',...
	  [0:20:70],'xlim',[2.8 7.4])
grid on
xl=xlabel('TriNet magnitude');
yl=ylabel('# of reporting stations');

% ERROR VERSUS NUMBER OF REPORTING STATIONS
axes(ah(3))
pav=plot(g.nonnans,abs(s-avhl),'o','Markerf','k','MarkerS',3,'MarkerE','k');
% Labels, etc
xl=xlabel('number of reporting stations');
yl=ylabel('|mag prediction error|');
grid on
% Cosmetics
hl=2;
set(ah(3),'ylim',[-0.3 1.3],'ytick',[-1:1],'xlim',[0 max(g.nonnans)+1])

fig2print(gcf,'portrait')
longticks(ah,2)

% Labels
[bh1,th1]=label(ah(1),'lr',12);
% Move this label about a bit
set(bh1,'verti',get(bh1,'verti')+[repmat(0.075,4,1) repmat(0,4,2)])
moveh(th1,0.075)
[bh2,th2]=label(ah(2:3),'ur',12,1);
set(bh2(1),'verti',get(bh2(1),'verti')+[repmat(0.075,4,1) repmat(0,4,2)])
set(bh2(2),'verti',get(bh2(2),'verti')+...
	   [repmat(1.06,4,1) repmat(0,4,2)])
moveh(th2(1),0.075)
moveh(th2(2),1.06)

figdisp

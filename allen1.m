function allen1
% ALLEN1
%
% Shows two seismograms and a zoom on the initial P-wave amplitude
% Watch that the panels sizes are exactly the same.
%
% Last modified by fjsimons-at-alum.mit.edu, 11/25/2008

% The commented material I really only needed once...

% Identify two earthquakes close of differing magnitude in some cluster
%[names,lon,lat,mag]=locident;
% Identify all the stations in another cluster 
%[statio,lon,lat]=statident([241.7 33.9],20);
% Identify two events of differing magnitude
% [mima,j]=min(mag);
% [mama,k]=max(mag);

% Set defaults
% defval('seis1',names(j,:))
% defval('seis2',names(k,:))
% disp(names(j,:))
% disp(names(k,:))

% Possible stations, from doing
% ls(fullfile(names(j,:),'BHZdata'))
% ls(fullfile(names(k,:),'BHZdata'))
% are:
% CWC, GSC, ISA
% stats=statident({'CWC','GSC','ISA'})
% Plot this, take the furthest, GSC

% Now I pretend I've done the above commented stuff

defval('ddir','/home/fjsimons/MyPapers/2006/EPSL-2006/DATA')
defval('seis1','20010730.233417')
defval('seis2','19950920.232736')
defval('mima',3.70)
defval('mama',5.76)

% Load data
[x1,h1,t1,p1,ts1]=readsac(fullfile(ddir,seis1,'BHZdata',...
				   'GSC.BHZ.sac.t0'),0);
[x2,h2,t2,p2,ts2]=readsac(fullfile(ddir,seis2,'BHZdata',...
				       'gsc.BHZ.sac.t0.sync'),0);

% Better check that the returned is what I think it is
difer(mima-h1.MAG,4)
difer(mama-h2.MAG,4)

% Figure out data amplification
mx1=range(minmax(x1));
mx2=range(minmax(x2));
xli1=[h1.T0-1 h1.T0+3];
sel1=ts1<xli1(2)&ts1>xli1(1);
mx3=range(minmax(x1(sel1)));
xli2=[h2.T0-1 h2.T0+3];
sel2=ts2<xli2(2)&ts2>xli2(1);
mx4=range(minmax(x2(sel2)));

% Make plot
clf
[ah,ha]=krijetem(subnum(3,2));
% Right edge of text in panel
txl=102;
% Ratio to scale white box underneath tex 
bs=0.65;

ylis=[-1.2 1.2];

leg1=sprintf('mag= %3.1f  az= %i%s%s= %3.1f%s',...
	     mima,round(h1.AZ),str2mat(176),...
	     '\Delta',h1.GCARC,str2mat(176));
leg2=sprintf('mag= %3.1f  az= %i%s%s= %3.1f%s',...
	     mama,round(h2.AZ),str2mat(176),...
	     '\Delta',h2.GCARC,str2mat(176));

axes(ah(1))
p(1)=plot(ts1,detrend(scale(x1),'constant'),'k');
xlim([0 100])
xl(1)=xlabel('time (s)');
yl(1)=ylabel(sprintf('amplitude %s %i','\times',round(mx2/mx1)));
hold on
pt(1)=plot([h1.T0 h1.T0],ylis,'k:');
% Watch for the quick fix here in the next line
ps(1)=plot([h1.T1 h1.T1],ylis+[0.2 0],'k:');
tx(1)=text(txl,-1.,leg1,'horizontala','right');
fb(1)=fillbox(ext2lrtb(tx(1),bs,bs),'w');
set(fb(1),'EdgeC','w')
top(tx(1),ah(1))

axes(ah(2))
p(1)=plot(ts2,detrend(scale(x2),'constant'),'k');
xlim([0 100])
xl(2)=xlabel('time (s)');
yl(2)=ylabel('amplitude');
hold on
pt(2)=plot([h2.T0 h2.T0],ylis,':');
ps(2)=plot([h2.T1 h2.T1],ylis,':');
tx(2)=text(txl,-1.,leg2,'horizontala','right');
fb(2)=fillbox(ext2lrtb(tx(2),bs,bs),'w');
set(fb(2),'EdgeC','w')
top(tx(2),ah(2))

axes(ah(3))
sec1=scale(x1(sel1));
p(1)=plot(ts1(sel1),sec1-sec1(1),'k');
xlim(xli1)
xl(3)=xlabel('time (s)');
yl(3)=ylabel(sprintf('amplitude %s %i','\times',round(mx4/mx3)));
hold on
pt(3)=plot([h1.T0 h1.T0],ylis,'k:');
px(1)=plot(ts1(sel1),zeros(size(ts1(sel1))),':');

axes(ah(4))
sec2=scale(x2(sel2));
p(1)=plot(ts2(sel2),sec2-sec2(1),'k');
xlim(xli2)
xl(4)=xlabel('time (s)');
yl(4)=ylabel('amplitude');
hold on
pt(4)=plot([h2.T0 h2.T0],ylis,'k:');
px(2)=plot(ts2(sel2),zeros(size(ts2(sel2))),'k:');

% Cosmetics
fig2print(gcf,'portrait')
delete(ah(5:6)); ah=ah(1:4); ha=ha([1 2 4 5]);
nolabels(ha,2); longticks(ah)
set(ah,'ylim',ylis)

% Somehow this messes it up as to the relative box sizes
%serre(ah(1:2),1/3,'across')
%serre(ah(3:4),1/3,'across')
%moveh(ah(1),.075)
%moveh(ah(3),.075)

% Note that dataaspectratio is different in both plots... so must give
% the posxmul etc explicitly
[a,b]=label(ah,[],12,[],[],[],[],0.75,0.75);

figdisp([],[],[],1)
% Remove the funky fake characters
!degs /home/fjsimons/EPS/allen1.eps

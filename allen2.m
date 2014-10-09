function allen2
% Time-frequency analysis of the seismograms studied in ALLEN1
%
% Last modified by fjsimons-at-alum.mit.edu, 11/30/2006

% Identify two earthquakes close of differing magnitude in some cluster
[names,lon,lat,mag]=locident;
% Identify all the stations in another cluster 
[statio,lon,lat]=statident([241.7 33.9],20);
% Identify two events of differing magnitude
[mima,j]=min(mag);
[mama,k]=max(mag);

% Set defaults
defval('ddir','/home/fjsimons/EALARMS')
defval('seis1',names(j,:))
defval('seis2',names(k,:))
disp(names(j,:))
disp(names(k,:))

% Load data
[x1,h1,t1,p1,ts1]=readsac(fullfile(ddir,seis1,'BHZdata',...
				   'GSC.BHZ.sac.t0'),0);
[x2,h2,t2,p2,ts2]=readsac(fullfile(ddir,seis2,'BHZdata',...
				       'gsc.BHZ.sac.t0.sync'),0);

% Lowpass filter seismograms at 10 Hz like Richard does
%x1f=lowpass(x1,1/h1.DELTA,10-1e-4);
%x2f=lowpass(x2,1/h2.DELTA,10-1e-4);
% Note that I checked this does not influence the waveform at all...
% in the P-wave window
% Just use Matlab's resample function, in that case
rrate=2;
x1f=resample(x1,1,rrate); 
ts1f=linspace(ts1(1),ts1(end),length(x1f));
x2f=resample(x2,1,rrate); 
ts2f=linspace(ts2(1),ts2(end),length(x2f));

% Replace with a synthetic:
% x1f=sin(2*pi*1.5*ts1f);
% x2f=sin(2*pi*2.5*ts2f);

% Note that I have checked this, too, it looks wonderful
% Make spectrogram of the first 100 seconds of the resampled data 
% That is, between 0 and 100!
sel1=ts1f<100&ts1f>0;
sel2=ts2f<100&ts2f>0;
% Sliding window of about this many s
wlens=4;
wl1=ceil(wlens./h1.DELTA/rrate);
wl2=ceil(wlens./h2.DELTA/rrate);
% Overlap, in percent
olap=50;
ol1=floor(olap/100*wl1);
ol2=floor(olap/100*wl2);
% Number of frequency sample
nfft=64;
[Ba2_1,F1,T1]=spectrogram(x1f(sel1),nfft,1./h1.DELTA/rrate,wl1,ol1);
[Ba2_2,F2,T2]=spectrogram(x2f(sel2),nfft,1./h2.DELTA/rrate,wl2,ol2);

% Now calculate the so-called instantaneous frequency
[in1,t1]=installen(x1f(sel1),ts1f(sel1),wlens,olap);
[in2,t2]=installen(x2f(sel2),ts2f(sel2),wlens,olap);
% The other measures are roughly similar... but noth worth showing. They
% are very noisy, and depend on a lot of smoothing... 

% Make plot
clf
[ah,ha]=krijetem(subnum(3,2));
colormap(flipud(gray(25)))

% Logarithmic amplitude scaling
los=[-30 0];

axes(ha(1))
imagesc(T1+wlens/2,F1,decibel(Ba2_1),los); axis xy
xlim([0 100])
hold on
pt(1)=plot([h1.T0 h1.T0],[-1 1]+minmax(F1),'k-');
ps(1)=plot([h1.T1 h1.T1],[-1 1]+minmax(F1),'k-');
% Plot the allen data at the midpoint of the windows
pa(1)=plot(t1+wlens/2,in1,'k','LineW',1);
yl(1)=ylabel('frequency (Hz)');
hold off

axes(ha(2))
imagesc(T2+wlens/2,F2,decibel(Ba2_2),los); axis xy
xlim([0 100])
colormap(flipud(gray(25)))
hold on
pt(2)=plot([h2.T0 h2.T0],[-1 1]+minmax(F1),'k-');
ps(2)=plot([h2.T1 h2.T1],[-1 1]+minmax(F2),'k-');
pa(2)=plot(t2+wlens/2,in2,'k','LineW',1);
xl(1)=xlabel('time (s)');
yl(2)=ylabel('frequency (Hz)');
hold off

cb=colorbarf('hor',10,'Helvetica',[0.1292  0.31  0.3367  0.0207]);
set(get(cb,'xlabel'),'string','spectral density (dB)')
longticks(cb,2)

% Cosmetics
fig2print(gcf,'portrait')
delete(ha(3:end))
ha=ha(1:2);
ah=ah([1 3]);
longticks(ah)
[a,b]=label(ah,'ur',12);
nolabels(ah(1),1)
serre(ha(1:2),1/3,'down')
set(ah,'ytick',0:5)

figdisp

disp('Mention in axis how much bigger the absolute values are')
disp('Compare dimensions with allen3')
disp(['Watch for data aspect ratio, correct by simply putting in' ...
      ' keyboard'])


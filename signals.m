function varargout=...
    signals(tipe,dec,ddir,nams,wsec,wolap,noxt,sax,fname,varargin)
% SIGNALS(tipe,dec,ddir,nams,wsec,wolap,noxt,sax,fname,fops)
% [Bl10,F,T,t,wlen,Fs,stt,idh]=SIGNALS(...)
%
% Makes diagnostic plots of (a pair of) "seismograms" 
% using SPECTRAL ANALYSIS:
%
% L  Time-domain representation of the signal
% M  Spectrogram: time-frequency representation
% R  Spectral density: using Chave's method
%
% If Header variables T0 and T1 are defined, plots those too.
%
% INPUT:
%
% tipe      0   uses up to two provided file name(s) in 'nams';
%          1-16 where signal/no signal plot pairs will be shown, of 
%          ----------------------------------------------------------
%           1 airguns                    2 blue whale call Type 1 
%           3 blue whale call Type 2     4 blue whale call Type 3
%           5 fin whale call             6 minke whale call
%           7 PKP-phase                  8 P-phase (local)
%           9 P,S,T phases (local)      10 P-phase (teleseismic)
%          11 ship noise                12 T-phase
%          13 P-phase                   14 P-phase
%          15 PKP-phase (also PKiKP)    16 PKP-phase (also PKiKP)
%          ----------------------------------------------------------
% dec      0 Nothing special happens to the data [default]
%          1 Reads decimated data from file *_dec.sac
%          2 Data in file will be filtered on the fly using defaults in 
%            the function TIMDOMPLOT
% ddir     Directory where data are located
% nams     Cell with up to TWO file name strings; if provided, tipe is
%          set to 0 - note this MUST be a cell array
% wsec     Window length in second
% wolap    Fraction of window overlap [0->1]
% noxt     0 Appends the .sac extension when looking for the filename [default]
%          1 Doesn't do this
% sax      1 Equalizes the y-axis limits in both rows of panels 
%          0 Doesn't do this [default]
% fname    Filter name (defaults taken from timdomplot)
% fops     Filter options (defaults taken from timdomplot)
%
% OUTPUT:
%
% Bl10     Spectral density in decibel from 10*log10(abs(B.^2))
%          Spectral density is energy per frequency (UNIT^2/HZ)
%          Thus, semilogx(flipud(sum(flipud(Bl10),2)/size(Bl10,2))) looks
%          like what specdensplot returns, except for an offset, since
%          you've been weighting different overlaps
% F        Frequency axis. The first AVAILABLE frequency F(1) is said to
%          be zero, but is really the Rayleigh frequency 1/N/Dt or 1 over
%          the window size indicated in the legend. The first VISIBLE
%          frequency F(2) on the logarithmic scale is the 0+fN/nfft*2,
%          i.e. the frequency discretization interval Fs/nfft. 
% T        Time axis of the spectrogram
% t        Axis handle of the title
% wlen     Number of samples
% Fs       Sampling frequency, in Hertz
% stt      Title handle
% idh      Plot id handle
% 
% See also: EVENTS, SIGNALS_ILL
% 
% EXAMPLE I
%
% signals(0,0,'/home/fjsimons/MERMAID/11042003/HRSECTIONS/REALDATA',...
%         {'mpilot1_001' 'mpilot1_002'})
%
% EXAMPLE II:
%
% for index=1:2:200
%   nams{1}=sprintf('mpilot3_%3.3i',index);
%   nams{2}=sprintf('mpilot3_%3.3i',index+1);
%   signals(0,0,'/home/fjsimons/MERMAID/08092007/HRSECTIONS/',nams)
% end
%
% EXAMPLE III:
% 
% signals(0,0,pwd,{'Okal_dec16'},6000)
%
% EXAMPLE IV:
%
% signals(3)
% 
% Last modified by fjsimons-at-alum.mit.edu, 09/13/2007

defval('tipe',1)
defval('dec',0)
defval('ddir','/home/fjsimons/MERMAID/SIGNALS');
defval('actpr',0) % Print, actually, in local directory
defval('wsec',5)
defval('wolap',0.875);
defval('noxt',0)
defval('sax',0)
defval('fname',[])

% I/O Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>=4
  tipe=0;
end
if tipe==0
  names{1}{1}=nams{1};
  if length(nams)==2
    names{1}{2}=nams{2};
    tits{1}=nounder(sprintf('%s / %s',nams{1},nams{2}));
  else
    tits{1}=nounder(nams{1});
  end
  tipe=1;
else
  nams=namit;
  names=namit;
  tits=titit;
end
if noxt==0
  sxt='.sac';
else
  sxt=[];
end
% End I/O Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
[ah,ha]=krijetem(subnum(2,3));
fig2print(gcf,'landscape')
idh=id; % Use fixcopy.pl later on!

% TIME-DOMAIN REPRESENTATION OF FIRST SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(1))
switch dec
 case {0,2}
  filenam1=fullfile(ddir,[names{tipe}{1} sxt]);
 case 1
  filenam1=fullfile(ddir,[names{tipe}{1} '_dec' sxt]);
end

[x1,h1,Fs,p(1),xl(1),yl(1),t(1)]=...
    timdomplot(filenam1,dec,fname,varargin{:});

% TIME-FREQUENCY SPECTROGRAM OF FIRST SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))
% The desired window length, in samples
wlen=floor(wsec*Fs);
% The number of frequencies, ideally the length of the window
nfft1=max(2^nextpow2(wlen),512);
[p(2),xl(2),yl(2),bm{1},Bl10,F,T]=...
    timspecplot(x1,h1,nfft1,Fs,wlen,wolap,h1.B);

% SPECTRAL DENSITY OF FIRST SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(3))
lwin=floor(h1.NPTS/2);
olap=70;
nfft2=512;
[p(3:6),xl(3),yl(3)]=specdensplot(x1,nfft2,Fs,lwin,olap);

% OUTPUT OF FIRST SIGNAL ONLY
if nargout>=6
  varargout{1}=Bl10;
  varargout{2}=F;
  varargout{3}=h1.B+wlen/Fs/2+T;
  varargout{4}=tits{tipe};
  varargout{5}=wlen;
  varargout{6}=Fs;
end

if prod(size(nams))>1
  % TIME-DOMAIN REPRESENTATION OF SECOND SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(4))
  switch dec
   case {0,2}
    filenam2=fullfile(ddir,[names{tipe}{2} sxt]);
   case 1
    filenam2=fullfile(ddir,[names{tipe}{2} '_dec' sxt]);
  end
  [x2,h2,Fs,p(7),xl(4),yl(4),t(2)]=timdomplot(filenam2,dec);

  % TIME-FREQUENCY SPECTROGRAM OF SECOND SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(5))
  [p(8),xl(5),yl(5),bm{2}]=timspecplot(x2,h2,nfft1,Fs,wlen,wolap,h2.B);

  % SPECTRAL DENSITY OF SECOND SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(6))
  [p(9:12),xl(6),yl(6)]=specdensplot(x2,nfft2,Fs,lwin,olap);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % COSMETICS FOR THE SECOND ROW OF PANELS
  if sax==1
    seemax(ha(1:2),2); seemax(ha(3:4),3); seemax(ha(5:6),2)
  end
  motion=-0.05;
  set(p([3 9]),'LineW',2,'Color','b')
  set(p([4 5 10 11]),'LineW',1,'Color',[0.7 0.7 0.7])
else
  motion=-0.2;
  set(p(3),'LineW',2,'Color','b')
  set(p([4 5]),'LineW',1,'Color',[0.7 0.7 0.7])
  delete(ah(4:6))
  ah=ah(1:3);
end  

% MORE COSMETICS APPLIED TO ALL PANELS
set(ah,'Xgrid','on','YGrid','on')
if ~isempty(t); delete(t); end
movev(ah(1:3),motion)
if ~isempty(bm{1}) & ~isempty(bm{2})
  bmstr1=repmat('%i ',1,length(bm{1}));
  bmstr2=repmat('%i ',1,length(bm{2}));
  stt=supertit(ah(1:3),sprintf(['%s \n' bmstr1 '/ ' bmstr2],...
			      tits{tipe},bm{1},bm{2}, ...
			      20));
else
  stt=supertit(ah(1:3),tits{tipe},20);
end
movev(stt,.3)
longticks(ah)

% PRINT OUTPUT SUGGESTIONS
switch dec
  case {0,2}
   if actpr==1
     print('-depsc',sprintf('signals_%s_%i',pref(names{1}{1}),dec))
   else
     % Set last flag to zero if you no longer want them printed
     figdisp(sprintf('signals_%s_%i',pref(names{1}{1}),dec),...
	     [],[],0)
     % Later, do fixcopy.pl on the files
   end
 case 1
  figdisp(sprintf('signals_%3.3i_dec',tipe))
end

if nargout>=7
  varargout{7}=stt;
end
if nargout>=8
  varargout{8}=idh;
end

% Return middle and range
%k=get(ah(1),'xlim');
%round(k(1)+range(k)/2)
%round(range(k)/2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nams=namit
nams{1}= {'airguns1'        'airguns2'};
nams{2}= {'bluewhale1_call' 'bluewhale1_nocall'};
nams{3}= {'bluewhale2_call' 'bluewhale2_nocall'};
nams{4}= {'bluewhale3_call' 'bluewhale3_nocall'};
nams{5}= {'finwhale_call'   'finwhale_nocall'};
nams{6}= {'minkewhale_call' 'minkewhale_nocall'};
nams{7}= {'PKP_call'        'PKP_nocall'};
nams{8}= {'Pphase_call'     'Pphase_nocall'};
nams{9}= {'PST_call'        'PST_nocall'};
nams{10}={'Ptele_call'      'Ptele_nocall'};
nams{11}={'shipnoise1'      'shipnoise2'};
nams{12}={'Tphase_call'     'Tphase_nocall'};
nams{13}={'P2_call'         'P2_nocall'};
nams{14}={'P3_call'         'P3_nocall'};
nams{15}={'PKiKP_call'      'PKiKP_nocall'};
nams{16}={'PKiKP2_call'     'PKiKP2_nocall'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tots=titit
tots{1}= {'Airgunning'};
tots{2}= {'Blue Whale calls I'};
tots{3}= {'Blue Whale calls II'};
tots{4}= {'Blue Whale calls III'};
tots{5}= {'Fin Whale calls I'};
tots{6}= {'Minke Whale calls I'};
tots{7}= {['PKP arrival (Mw 7.4 at 153.0' sprintf('%s',str2mat(176)) ')']};
tots{8}= {['Local P arrival (at 6.2' sprintf('%s',str2mat(176)) ')']};
tots{9}= {['P, S and T phases (Mw 5.5 at 3.6' sprintf('%s',str2mat(176)) ')']};
tots{10}={['Teleseismic P arrival (Mw 7.0 at 82.9' sprintf('%s',str2mat(176)) ')']};
tots{11}={'Ship Noise'};
tots{12}={['T phase (Mw 5.5 at 10.9' sprintf('%s',str2mat(176)) ')']};
tots{13}={['Teleseismic P arrival (Mw 6.9 at 66.0' sprintf('%s',str2mat(176)) ')']};
tots{14}={['Teleseismic P arrival (Mw 6.9 at 55.7' sprintf('%s',str2mat(176)) ')']};
tots{15}={['PKP arrival (Mw 7.9 at 145.8' sprintf('%s',str2mat(176)) ')']};
tots{16}={['PKP arrival (Mw 7.9 at 131.5' sprintf('%s',str2mat(176)) ')']};

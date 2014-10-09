function varargout=signals2(tipe,dec,ddir,nams,wsec,wolap,noxt,sax,...
			    tipo,nvm,nd,pph,intel,thresh)
% [Bl10,F,T,t,wlen,Fs,stt,idh]=SIGNALS(tipe,dec,ddir,nams,wsec,wolap,noxt,sax,...
%			    tipo,nvm,nd,pph,intel,thresh)
%
% Makes a plot that relies on combining WTEVENTS with SIGNALS
% 
% INPUT:    
% 
% tipe,dec,ddir,nams,wsec,wolap,noxt,sax     As listed in SIGNALS
% nvm,nd,pph,intel,thresh                    As listed in WTEVENTS
%
% OUTPUT:
%
% Bl10,F,T,t,wlen,Fs,stt,idh                 As listed in SIGNALS
%
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2010

defval('tipe',1)
defval('dec',0)
defval('ddir','/home/fjsimons/MERMAID/SIGNALS');
defval('actpr',0) % Print, actually
defval('wsec',5)
defval('wolap',0.875);
defval('noxt',0)
defval('sax',1)
defval('tipo','CDF')
defval('nvm',[2 4])
defval('nd',5)
defval('pph',3)
defval('intel',0)
defval('thresh',0)

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
    tits{1}=nams{1};
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

% TIME-DOMAIN REPRESENTATION OF FIRST SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ha(1))
switch dec
 case {0,2}
  filenam1=fullfile(ddir,[names{tipe}{1} sxt]);
 case 1
  filenam1=fullfile(ddir,[names{tipe}{1} '_dec' sxt]);
end
[x1,h1,Fs,p(1),xl(1),yl(1),t(1)]=timdomplot(filenam1,dec);

% TIME-FREQUENCY SPECTROGRAM OF FIRST SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))
% The desired window length, in samples
wlen=floor(wsec*Fs);
% The number of frequencies, ideally the length of the window
nfft1=max(2^nextpow2(wlen),512);
[p(2),xl(2),yl(2),bm{1},Bl10,F,T]=...
    timspecplot(x1,h1,nfft1,Fs,wlen,wolap,h1.B);

% WAVELET ANALYSIS OF FIRST SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(3))
[psg{1},xl(4),yl(4)]=scalogramplot(x1,h1,tipo,nvm,nd,pph,intel,thresh);

% OUTPUT OF FIRST SIGNAL ONLY
if nargout==6
  varargout{1}=Bl10
  varargout{2}=F;
  varargout{3}=h1.B+wlen/Fs/2+T;
  varargout{4}=tits{tipe};
  varargout{5}=wlen;
  varargout{6}=Fs;
end

if prod(size(nams))>1
  % TIME-DOMAIN REPRESENTATION OF SECOND SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ha(2))
  switch dec
   case {0,2}
    filenam2=fullfile(ddir,[nams{tipe}{2} '.sac']);
   case 1
    filenam2=fullfile(ddir,[nams{tipe}{2} '_dec.sac']);
  end
  [x2,h2,Fs,p(7),xl(4),yl(4),t(2)]=timdomplot(filenam2,dec);
  
  % TIME-FREQUENCY SPECTROGRAM OF SECOND SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(5))
  [p(8),xl(5),yl(5),bm{2}]=timspecplot(x2,h2,nfft1,Fs,wlen,wolap,h2.B);
  
  % WAVELET ANALYSIS OF SECOND SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(6))
  [psg{2},xl(8),yl(8)]=scalogramplot(x2,h2,tipo,nvm,nd,pph,intel,thresh);

  % COSMETICS FOR THE SECOND ROW OF PANELS
  if sax==1
    seemax(ha(1:2),2); seemax(ha(3:4),3) ; seemax(ha(5:6),3)
  end
  motion=-0.05;
else
  motion=-0.2;
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
  st=supertit(ah(1:4),sprintf(['%s \n' bmstr1 '/ ' bmstr2],...
			      tits{tipe},bm{1},bm{2}, ...
			      20));
else
  st=supertit(ah(1:3),tits{tipe},20);
end
movev(st,.3)
%moveh(ha(5:6),-.0125)
longticks(ah)

% PRINT OUT SUGGESTIONS
% Need to do -zbuffer if not Postscript chokes on imagesc "rangecheck"
switch dec
  case {0,2}
   if actpr==1
     print('-depsc','-zbuffer','-r600',...
	   sprintf('/home/fjsimons/GifPix/EPS/signals2_%3.3i',tipe))
   else
     figdisp(sprintf('signals2_%3.3i',tipe))
   end
 case 1
  if actpr==1
    print('-depsc','-zbuffer','-r600',...
	  sprintf('/home/fjsimons/GifPix/EPS/signals2_%3.3i_dec',tipe))
  else
    figdisp(sprintf('signals2_%3.3i_dec',tipe))
  end
end

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
tots{7}= {['PKP arrival (Mw 7.4 at 153.0' str2mat(176) ')']};
tots{8}= {['Local P arrival (at 6.2' str2mat(176) ')']};
tots{9}= {['P, S and T phases (Mw 5.5 at 3.6' str2mat(176) ')']};
tots{10}={['Teleseismic P arrival (Mw 7.0 at 82.9' str2mat(176) ')']};
tots{11}={'Ship Noise'};
tots{12}={['T phase (Mw 5.5 at 10.9' str2mat(176) ')']};
tots{13}={['Teleseismic P arrival (Mw 6.9 at 66.0' str2mat(176) ')']};
tots{14}={['Teleseismic P arrival (Mw 6.9 at 55.7' str2mat(176) ')']};
tots{15}={['PKP arrival (Mw 7.9 at 145.8' str2mat(176) ')']};
tots{16}={['PKP arrival (Mw 7.9 at 131.5' str2mat(176) ')']};



function varargout=wtevents3(findex,ddir,tipe,nvm,nd,pph,intel)
% [qtf,qts,qualts,qualtf]=WTEVENTS3(findex,ddir,tipe,nvm,nd,pph,intel)
%
% Compares detection quality analysis from the time-frequency (BFT) and the
% time-scale (BST) methods.
%
% INPUT:
%
% findex           [scalar] a running number identifying the file in ddir
%                  [string] a filename string
% ddir             is the directory (with the events identified by STALTA)
%                  [default:~/MERMAID/11042003/EVENTS]
% tipe             wavelet type ['CDF']
% nvm              number of vanishing moments [[2 4]]
% nd               number of cascades
% pph              1 Time-domain full rate
%                  2 Time-domain polyphase
%                  3 Z-domain polyphase
%                  4 Lifting implementation
% intel            1 With integer rounding (for lifting only)
%                  0 Without integer rounding
%
% SEE ALSO: EVENTS, SIGNALS, WTEVENTS, WTEVENTS2, BFT, BST
%
% EXAMPLE:
%
% wtevents3('demo1') % Runs through the positive MERMAID detections
% wtevents3('demo2') % Runs through some other ones
%
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2010

% Two nice examples are: for 3 and 4

% Establish defaults for the data location
defval('findex',1)
defval('ddir','/home/fjsimons/MERMAID/09102004/EVENTS/')

if ~isstr(findex) || isempty(strmatch('demo',findex))
  % Establish defaults for the wavelet analysis
  defval('tipe','CDF')
  defval('nvm',[2 4])
  defval('nd',5)
  defval('pph',3)
  defval('intel',0)

  tph{1}='Time-Domain Full Rate';
  tph{2}='Time-Domain Polyphase';
  tph{3}='Z-Domain Polyphase';
  tph{4}='Lifting Implementation';

  tph=tph{pph};

  % Get data
  if ~isstr(findex)
    files=ls2cell(ddir);
    filen=files{findex};
  else
    filen=findex;
  end

  clf
  ah=krijetem(subnum(2,3)); 
  delete(ah(4:6)); ah=ah(1:3);

  % FILTERED TIME-DOMAIN REPRESENTATION OF FIRST SIGNAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  filenam1=fullfile(ddir,filen);
  [x1,h1,t{1},p(1)]=readsac(filenam1,0);
  cof=2;
  Fs=1/h1.DELTA;
  sig=lowpass(x1,Fs,cof,2,2,'butter');
  timax=linspace(h1.B,h1.E,h1.NPTS);
  p(1)=plot(timax,sig);
  xll(1)=xlabel(sprintf('%s ; %i s selected',...
			'Time (s)',ceil(h1.T1-h1.T0)));
  axis tight
  yl(1)=ylabel(sprintf('Filtered Amplitude, %i Hz',cof));
  hold on
  grid on
  yli=ylim;
  pt01(:,1)=plot(repmat([h1.T0 h1.T1],2,1),[yli(1) yli(2)],'k--');

  % TIME-FREQUENCY ANALYSIS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  wolap=0.875;
  % Window length needs to be smaller than nfft
  nfft2=512;
  wlen=min(nfft2,floor(5*Fs));
  [B,F,T]=spectrogram(x1,nfft2,Fs,wlen,ceil(wolap*wlen));
  disp(sprintf('Window size for spectrogram:    %8.1f s',wlen/Fs))
  p(2)=imagesc(h1.B+wlen/Fs/2+T,F,10*log10(B));
  hold on
  yli=ylim;
  pt01(:,2)=plot(repmat([h1.T0 h1.T1],2,1),[yli(1) yli(2)],'k--');
  axis xy; colormap(jet)    

  % TIME-SCALE ANALYSIS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(3))
  ochk=1*mod(length(x1),2);
  [a,d,an,dn,xs]=succapp(tipe,nvm,nd,x1(1:end-ochk),pph,intel);
  [tsamp,tskol]=wtimax(a,d,an,length(x1),[h1.B h1.E]);
  meth=1;
  col='bw';
  [pd{3},stp,rgp]=dyadplot(x1(1:end-ochk),a,d,an,dn,meth,[h1.B h1.E],col,3);
  hold on
  yli=ylim;
  pt01(:,2)=plot(repmat([h1.T0 h1.T1],2,1),[yli(1) yli(2)],'k--');
  if meth==1 & col ~= 'bw'
    colormap(flipud(gray(128)))
    caxis([0 3*stp])
  end
  set(ah(3),'Xlim',[h1.B h1.E])

  % QUALITY ANALYSIS
  axes(ah(2))
  [qualtf,qtf]=bft(10*log10(B),F,h1.B+wlen/Fs/2+T,h1.T0,h1.T1);
  xl(2)=xlabel(sprintf('%s ; %3.1f s window, SN %5.1f','Time (s)',...
		       wlen/Fs,qtf));
  yl(5)=ylabel(sprintf('%s ; Nyquist %5.1f Hz','Frequency (Hz)',Fs/2));
  tl(1)=title(sprintf('%s %i','Quality',qualtf));

  axes(ah(3))
  % Feed it wavelet coefficients and times
  [qualts,qts]=bst(d,tskol,h1.T0,h1.T1);
  yll(3)=ylabel('Resolution level (Coarseness)');
  xll(3)=xlabel(sprintf('%s ; SN %5.1f','Time (s)',qts));
  tl(2)=title(sprintf('%s %i','Quality',qualts));

  % And judge the results
  titl=sprintf('%s / %s %i,%i / %s',filen,tipe,nvm,tph);
  titl(find(abs(titl)==95))='-';

  st=supertit(ah,titl);
  movev(ah,-.2)
  movev(st,-1.75)

  if pph==3
    if isstr(findex)
      figdisp(sprintf('wtevents3_%s',pref(findex)))
    else
      figdisp(sprintf('wtevents3_%3.3i',findex))
    end
  elseif pph==4
    if isstr(findex)
      figdisp(sprintf('wtevents3_%s_l',pref(findex)))
    else
      figdisp(sprintf('wtevents3_%3.3i_l',findex))
    end
  end
  
  % Output
  varns={qtf,qts,qualts,qualtf};
  varargout=varns(1:nargout);
elseif strcmp(findex,'demo1')
  % These are the positive identifications from our experiments and
  % should all get a quality criterion of 1 in the wavelet space
  dirs={'11042003','11042003','11042003',...
	'09102004','08092007'};
  files={'mpilot1_21118','mpilot1_98947','mpilot1_138106',...
	 'mpilot2_79891','mpilot3_4759'};
  for ind=1:length(dirs)
    wtevents3(sprintf('%s.sac',files{ind}),...
	      sprintf('/home/fjsimons/MERMAID/%s/EVENTS/',dirs{ind}))
    keyboard
  end
elseif strcmp(findex,'demo2')
  % These are some sections from other hydrophones
  % While the code does not break, we certainly shouldn't run it as is
  % without filtering, etc. as these data sections are very high frequency
  files={'16n043w_1','16n043w_2','16n049w_1', ...
	 '16n049w_2','26n040w_1','26n040w_2', ...
	 '26n050w_1','26n050w_2','32n035w_1', ...
	 '32n035w_2','35n043w_1','35n043w_2 '};
  for ind=1:length(files)
    wtevents3(sprintf('%s.sac',files{ind}),...
	      sprintf('/home/fjsimons/MERMAID/SIGNALS/'))
    keyboard
  end
end


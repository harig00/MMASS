function varargout=bayram(X,O,nfft,sint,R,L,X2)
% BAYRAM(X,O,nfft,sint,R,L)
% [tim,freq,SX,H,V]=BAYRAM(X,O,nfft,sint,R,L)
% [tim,freq,SX,H,V,SX2,SX12,C2]=BAYRAM(X1,O,nfft,sint,R,L,X2)
%
% This program takes a time series and calculates the Hermite-multiple
% window spectral estimate of it (including coherence for two signals)
%
% X      analyzed signal (or X1 and X2)
% O      overlap of the windows in percent
% nfft   frequency sampling for each of the ffts of the windows
% sint   time sampling interval DX (default is 1)
% R      radius of conentration in the phase plane
% L      time resolution of the analysis (window length in samples)
% X2     second signal in case coherence is required
%
% tim    time domain  of the signal
% freq   frequency domain
% SX     spectral estimate (or SX and SY and C2)
% H      Hermite matrix used
% V      matrix with eigenvalues
%
% See also COHEN
%
% EXAMPLE 1
%
% nfft=1024; t=linspace(0,1,nfft); sint=indeks(diff(t),1);
% X=chirp(t,0,1,200,'quadratic');
% X=sigmerge(X,rand(size(X)),20);
% subplot(211); plot(t,X)
% subplot(212);
% bayram(X,99,256,sint,2);
%
% EXAMPLE 2
%
% X=rand(2048,1); Y=rand(2048,1); sint=1; t=1:2048; R=3; L=128;
% [tim,freq,SX,H,V,SX2,SX12,C2]=bayram(X,99,256,sint,R,L,Y);
% imagesc(tim,freq,C2); axis xy; caxis([0 1]);
%
% EXAMPLE 3
%
% subplot(211);  [y,t,sint]=notch; subplot(212)
% bayram(y,99,256,sint,2)

% filename= 'lowry1';
% dapa= '/rosella/data14/fjsimons/GRAVITY/READYSETS';
% load(fullfile(dapa,filename));

% Bayram and Baraniuk, Proc. IEEE Sig. Proc. 173-176, 1996

% One or two signals?
if  exist('X2') & ~isempty(X2)
  flag=1; else; flag=0;
end

% Preliminaries - demean signals
X=X(:)-mean(X); if flag; X2=X2(:)-mean(X2); end
if flag ; if length(X)~=length(X2); error('Data vectors do not match'); end;
end
nx=length(X);
if isempty(nfft); nfft=nx; end
if isempty(sint); sint=1; end
[H,V,K]=choice(R,L);

disp([ 'Number of Hermite windows: ',num2str(K),'; Time window ',sprintf('%4.2f',L*sint)])

[fax,select]=fftaxis1D(X,nfft,sint*(nx-1));
O=floor(O/100*L);
NWI=fix((nx-O)/(L-O)); % Number of windows
disp([ 'Number of time windows: ',num2str(NWI)])
colindex = 1 + (0:(NWI-1))*(L-O);
tim = (colindex-1)'*sint+(L*sint/2); % Middle time of window
rowindex = (1:L)';
if length(X)<(L+colindex(NWI)-1)
  keyboard
  X(L+colindex(NWI)-1) = 0;   % zero-pad last section of X
  if flag ; X2(L+colindex(NWI)-1) = 0; end
end
Ypre= zeros(L,NWI);  if flag; Y2pre=Ypre; end
% Put X into columns of Y with the proper offset
Ypre= X(rowindex(:,ones(1,NWI))+colindex(ones(L,1),:)-1);
if flag ; Y2pre= X2(rowindex(:,ones(1,NWI))+colindex(ones(L,1),:)-1); end

% Apply the windows to the array of offset signal segments.
SX=zeros(length(select),NWI); if flag; SX2=SX; SX12=SX;end
for index=1:K
  WI=H(:,index);
  Y = WI(:,ones(1,NWI)).*Ypre;
  if flag; Y2 = WI(:,ones(1,NWI)).*Y2pre; end

  % Now FFT Y; this does the columns
  FY = fft(Y,nfft); if flag; FY2= fft(Y2,nfft); end
  if isreal(X) % Assume both signals are both REAL or COMPLEX
    FY = FY(select,:); if flag; FY2=FY2(select,:); end
  end
  % Make weighted sum
  SX=V(index)*(abs(FY).^2)+SX; 
  if flag; SX2=V(index)*(abs(FY2).^2)+SX2; end
  if flag; SX12=V(index)*(FY.*conj(FY2))+SX12; end
end
SX=SX/sum(V); 
if flag
  SX2=SX2/sum(V); 
  SX12=SX12/sum(V); 
  C2=abs(SX12).^2./SX./SX2;
end

% Which property is plotted?
prop=SX;
%prop=decibel(SX);

if ~nargout
  axes(gca)
  imagesc(tim,fax,prop);
  axis xy; colormap(jet)  
  xlim([0 (nx-1)*sint]);
  xlabel('Time')
  ylabel('Frequency')
  title([ 'Number of Hermite windows: ',num2str(K),'; Time window ',sprintf('%4.2f',L*sint)])
else
  posout={ 'tim','fax','SX','H','V','SX2','SX12','C2'};
  for index=1:nargout
    eval([ 'varargout{index}=',posout{index},';'])
  end
end


function [FX,FY,SX,SY,SXY,COH2,E,W,covar,ADM]=mtm3(X,Y,NW,k)
% [FX,FY,SX,SY,SXY,COH2,E,W,covar,ADM]=MTM3(X,Y,NW,k)
%
% This is a memory efficient version of the retired two-dimensional
% multi-taper spectral analysis program MTM2. Speed is compromised.
% Made adaptation to 1-D time series; returns estimation variance
% Low wavenumbers (long wavelengths, dc-term) in the center
% If 'k' is not defined then takes k=max(round(2*NW)-1,1), 
% i. e. the Shannon number-1
%
% The power density spectra are SX, SY, and SXY.
%
% SEE ALSO: MTM
%
% Last modified by fjsimons-alum-mit.edu, 11/28/2010

disp('MTM3')
defval('NW',3)
defval('k',max(round(2*NW)-1,1)); % Shannon number-1

% Error check
if ~all(size(X)==size(Y)) 
  if prod(size(X).*size(Y))
    error('If both matrices are not empty, must be equal in size')
  end
end

% In case time series, keep only one dimension - make it a column vector
if ~any(size(X)==[1 1])
  % Detrend data with function 'planefit'
  disp('Data sets detrended')
  [a,b,c,d,e,Xfit]=planefit(X);
  [a,b,c,d,e,Yfit]=planefit(Y);
  X=X-Xfit;
  Y=Y-Yfit;
  fafa=0;
else
  disp('It''s a row !')
  X=detrend(X(:)');
  Y=detrend(Y(:)');
  fafa=1;
end
  
% Demean data
disp('Data sets demeaned')
X=X-nanmean(X(:));
Y=Y-nanmean(Y(:));

[irow,icol]=size(X);
disp([ 'Duration x half bandwidth product = ',num2str(NW)])
disp([ 'Shannon number = ',num2str(2*NW)])
disp([ 'Number of tapers used in 1D = ',num2str(k)])

kor=k;

% Resolution of the FFT-routine
nfftrow=irow;
nfftcol=icol;

% Compute the data windows for default parameters - 
% the windows are the columns of E
if fafa==0
  [Erow,Vrow]=dpss(irow,NW,k);
  Erow=Erow(:,1:k); Vrow=Vrow(1:k);
end
[Ecol,Vcol]=dpss(icol,NW,k);
Ecol=Ecol(:,1:k); Vcol=Vcol(1:k);

if fafa==0
  % Need all combinations of row and column windows.
  dims=ndims(X);
  [i j]=ind2sub([k k],1:k^2); 
  indises=[i ; j]';
  weights=Vcol(indises);
  k=length(indises);
end 

% Initializing arrays makes it sooo much faster
fX=repmat(NaN,[nfftrow nfftcol k]);
if ~isempty(Y); fY=repmat(NaN,[nfftrow nfftcol k]); end

plots=0;
if plots==1
  clf
  % This is the tapered data
  surf(X); shading flat; view(45,0)
end

if fafa==0
  % Create an 3-dimensional matrix ijk with k=1,...,k as many as there are
  % 2-dimensional data windows.  Calculate 2-dimensional FFT
  for index=1:k  
    fX(:,:,index)=fft(...
	repmat(Erow(1:irow,indises(index,1)),1,nfftcol).*...
	fft(repmat(Ecol(1:icol,indises(index,2)),1,irow)'.*X,nfftcol,2),...
	nfftrow,1);
    if ~isempty(Y)
      fY(:,:,index)=fft(...
	  repmat(Erow(1:irow,indises(index,1)),1,nfftcol).*...
	  fft(repmat(Ecol(1:icol,indises(index,2)),1,irow)'.*Y,nfftcol,2),...
	  nfftrow,1);
    end
  end
  
  % Calculate Fourier Transform--------------------------------------------
  % Sum everything up along third dimension in weighted fashion
  for index=1:k
    FX(:,:,index)=fX(:,:,index)*prod(weights(index,1:2));
    if ~isempty(Y) ; FY(:,:,index)=fY(:,:,index)*prod(weights(index,1:2)); end 
  end
  FX=fftshift(sum(FX,3)/sum(prod(weights,2)));
  if ~isempty(Y) ; FY=fftshift(sum(FY,3)/sum(prod(weights,2))); end 
  
  % Calculate Power Spectral Density--------------------------------
  for index=1:k
    SX(:,:,index)=(fX(:,:,index).*conj(fX(:,:,index)))...
	*prod(weights(index,1:2));
    if ~isempty(Y) ; SY(:,:,index)=(fY(:,:,index).*conj(fY(:,:,index)))...
	*prod(weights(index,1:2)); 
    end 
  end
  SX=fftshift(sum(SX,3)/sum(prod(weights,2)));
  if ~isempty(Y) ; SY=fftshift(sum(SY,3)/sum(prod(weights,2))); end 
  
  % Calculate Cross Spectral Density---------------------------------
  if ~isempty(Y)
    for index=1:k
      SXY(:,:,index)=(fX(:,:,index).*conj(fY(:,:,index)))...
	  *prod(weights(index,1:2));
    end
    SXY=fftshift(sum(SXY,3)/sum(prod(weights,2)));
  end
else
  % For one-dimensional signals
  % Only half of the returns will be meaningful
  % Include DC AND Nyquist
  select=[1:floor(nfftcol/2)+1];

  Xwigs=fft(Ecol.*repmat(X,k,1)',nfftcol);

  FX=Xwigs*Vcol;
  SX=(Xwigs.*conj(Xwigs))*Vcol;

  FX=FX(select);
  SX=SX(select);
  
  if ~isempty(Y)
    Ywigs=fft(Ecol.*repmat(Y,k,1)',nfftcol);
    FY=Ywigs*Vcol;
    SY=(Ywigs.*conj(Ywigs))*Vcol;
    SXY=(Xwigs.*conj(Ywigs))*Vcol;

    FY=FY(select);
    SY=SY(select);
    SXY=SXY(select);
  end
end
  
% Calculate Coherence---------------------------------------------
if ~isempty(Y)
  COH2=abs(SXY).^2./SX./SY;
  ADM=SXY./SX;
  % Calculate variance of coherence, Thomson, 1980, p 1087
  % REPLACE THIS WITH THE CRAMER-RAO BOUND
  % FOR THE VARIANCE OF THE COHERENCE SQUARED, WHICH IS NOT A NORMALLY
  % DISTRIBUTED VARIABLE...., THIS IS WHY WE DON'T SIMPLY TAKE THE SAMPLE
  % VARIANCE. But this is a lower bound...
  if fafa==0
    % covar=(k-1)/k/(k+1)*(1/k+2*(k-2)/(k+2).*sqrt(COH2));
    covar=2*COH2.*(1-COH2).^2/k;
  else
    % covar=(kor-1)/kor/(kor+1)*(1/kor+2*(kor-2)/(kor+2).*sqrt(COH2));
    covar=2*COH2.*(1-COH2).^2/kor;
  end
end

if isempty(Y) ; FY=[] ; SY=[] ; COH2=[]; SXY=[] ; covar=[] ; end

E=Ecol;

if fafa==0
  W=prod(weights,2);
else
  W=Vcol;
end

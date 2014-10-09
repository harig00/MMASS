function [FX,FY,SX,SY,SXY,COH2,E,W,coh2var,ADM]=mtm4(X,Y,R)
% [FX,FY,SX,SY,SXY,COH2,E,W,coh2var,ADM]=MTM4(X,Y,R)
% [FX,FY,SX,SY,SXY,COH2,E,W,coh2var,ADM]=MTM4(X,Y,R)
%
% Like MTM3 but with Hermite function windows for 2-D.
% 'R' is the concentration radius in the time domain, in seconds.
% Becaus the support goes form -5 to 5 in CHOICE, R/10 gives
% the fractional concentration in the time window. R=3 works well.
% Overlap can be adjusted to accomodate energy loss outside abs(t)>R.
% Admittance is what we need if the second field is gravity.
%
% Weighting is with the EIGENVALUE, WHICH IS GOOD FOR 
% WHITE SPECTRA (PW 369) BUT ADAPTIVE WEIGHTINGS EXIST

% Last modified by fjsimons-at--alum-mit.edu, Dec 18th, 2001

defval('R',3)

disp('MTM4 - Hermite method')

% Error check
if ~all(size(X)==size(Y)) 
  if prod(size(X).*size(Y))
    error('If both matrices are not empty, must be equal in size')
  end
end

% Detrend data with function 'planefit'
if ~any(size(X)==[1 1])
  % Detrend data with function 'planefit'
  disp('Data sets detrended')
  X=X-planefit(X);
  if ~isempty(Y) ; Y=Y-planefit(Y) ; end
  fafa=0;
else
  disp('It''s a row !')
  X=detrend(X(:)');
  if ~isempty(Y) ; Y=detrend(Y(:)') ; end
  fafa=1;
end
  
% Demean data
disp('Data sets demeaned')
X=X-nanmean(X(:));
if ~isempty(Y)
  Y=Y-nanmean(Y(:)) 
end

[irow,icol]=size(X);

% Choose Hermite windows
if fafa==0
  [Erow,Vrow,k,t]=choice(R,irow);
end
[Ecol,Vcol,k,t]=choice(R,icol);

% To illustrate the need
%Erow=ones(size(Erow));
%Ecol=ones(size(Ecol));

% For comparison with Slepian multitaper 
% normalize with sampling interval 
Erow=Erow*sqrt(10/irow);
Ecol=Ecol*sqrt(10/icol);

disp(sprintf(...
    [ 'Normalization:' repmat('%4.3f ',1,ceil(R^2))],...
    diag(Erow'*Erow)))
disp([ 'Concentration radius = '  ,num2str(R)])
disp([ 'Number of tapers used in 1D = ',num2str(k)])

kor=k;

% Resolution of the FFT-routine
nfftrow=irow;
nfftcol=icol;

% Need all combinations of row and column windows.
if fafa==0
  dims=ndims(X);
  [i j]=ind2sub([k k],1:k^2); 
  indises=[i ; j]';
  weights=Vcol(indises);
  k=size(indises,1);
  disp([ 'Number of tapers used in 2D = ',num2str(k)])
end 
  
% Initializing arrays makes it sooo much faster.
fX=repmat(NaN,[nfftrow nfftcol k]); FX=fX; SX=fX;
if ~isempty(Y)
  fY=repmat(NaN,[nfftrow nfftcol k]); FY=fY; SY=fY; SXY=fY
end

if fafa==0
  % Create an 3-dimensional matrix ijk with k=1,...,k as many as there are
  % 2-dimensional data windows.  Calculate 2-dimensional FFT
  for index=1:k  
    %  disp(sprintf('Calculating 2-D fft nr. %2d',index))
    fX(:,:,index)=fft(...
	repmat(Erow(1:irow,indises(index,1)),1,nfftcol).*...
	fft(repmat(Ecol(1:icol,indises(index,2)),1,irow)'.*...
	    X,nfftcol,2),nfftrow,1);
    if ~isempty(Y)
      fY(:,:,index)=fft(...
	  repmat(Erow(1:irow,indises(index,1)),1,nfftcol).*...
	  fft(repmat(Ecol(1:icol,indises(index,2)),1,irow)'.*...
	      Y,nfftcol,2),nfftrow,1);
    end
  end
  
  % Calculate Fourier Transform--------------------------------------------
  % Sum everything up along third dimension in weighted fashion
  for index=1:k
    FX(:,:,index)=fX(:,:,index)*prod(weights(index,1:2));
    if ~isempty(Y) 
      FY(:,:,index)=fY(:,:,index)*prod(weights(index,1:2));
    end 
  end
  FX=fftshift(sum(FX,3)/sum(prod(weights,2)));
  if ~isempty(Y) 
    FY=fftshift(sum(FY,3)/sum(prod(weights,2)));
  end 
  
  % Calculate Power Spectral Density--------------------------------
  for index=1:k
    SX(:,:,index)=(fX(:,:,index).*conj(fX(:,:,index)))...
	*prod(weights(index,1:2));
    if ~isempty(Y) ; SY(:,:,index)=(fY(:,:,index).*conj(fY(:,:,index)))...
	  *prod(weights(index,1:2)); 
    end 
  end
  SX=fftshift(sum(SX,3)/sum(prod(weights,2)));
  if ~isempty(Y) 
    SY=fftshift(sum(SY,3)/sum(prod(weights,2))); 
  end 
    
  % Calculate Cross Spectral Density---------------------------------
  if ~isempty(Y)
    for index=1:k
      SXY(:,:,index)=(fX(:,:,index).*conj(fY(:,:,index)))...
	  *prod(weights(index,1:2));
    end
    SXY=fftshift(sum(SXY,3)/sum(prod(weights,2)));
  end
else
  % include DC AND Nyquist
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
    % FOR THE VARIANCE OF THE COHERENCE SQUARED
    if fafa==0
      % coh2var=(k-1)/k/(k+1)*(1/k+2*(k-2)/(k+2).*sqrt(COH2));
      coh2var=2*COH2.*(1-COH2).^2/k;
    else
      % coh2var=(kor-1)/kor/(kor+1)*(1/kor+2*(kor-2)/(kor+2).*sqrt(COH2));
      coh2var=2*COH2.*(1-COH2).^2/kor;
    end
end

if isempty(Y) ; FY=[] ; SY=[] ; COH2=[] ; SXY=[] ; coh2var=[] ; end

if fafa==0
  E=squeeze(Ecol(:,1,:));
  W=weights;
else
  E=Ecol;
  W=Vcol;
end


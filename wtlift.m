function [a,d,xout,err]=wtlift(x,tipe,N,M,intel)
% [a,d,xout,err]=WTLIFT(x,tipe,N,M,intel)
%
% Fast discrete wavelet transform by (integer) lifting.
% SINGLE iteration of the filter bank!
%
% INPUT:
%
% x       Input signal
% tipe    Type of wavelet used ('CDF')
% N       Number of vanishing moments (primal step)
% M       Number of vanishing moments (dual step)
% intel   1 for INTEGER lifting
%         0 for REGULAR lifting
%
% OUTPUT:
%
% a       Cell array with scaling coefficients
% d       Cell array with wavelet coefficients
% xout    Reconstructed signal after forward/inverse transform pair
% err     RMSE of the reconstructed signal
%
% EXAMPLE:
%
% x=real(dopnoise(1024,200,60,10,70,128));
% [a,d,xout,err]=wtlift(x+rand(size(x))+1);
%
% See also LIFTCO, FBANK, ABANK, SBANK.

t0=clock;

defval('intel',1)
defval('tipe','CDF')
defval('N',2)
defval('M',4)
defval('x',[0 0 0 0 1 -1 2 5 1 7 0 0 0 0 0 0]')

x=x(:);

% Find wavelet coefficients in the catalog
[h0,f0,Pa,Ua,Kp,Ku]=wc(tipe,[N M]);
% In the case of CDF [2,4] we get:
% Pa=[1/2 1/2];
% Ua=[-3 19 19 -3]/64;
% Kp=sqrt(2)
% Ku=1/sqrt(2);
% and we need no more from WC.

% If using integer lifting, must have scale factors as lifting steps
if intel==1 & ~strcmp(tipe,'CDFI')
  disp('No precomputed lifting steps stored; attempt to reconstruct') 
  % Here we now need to absorb the scaling as four more lifting steps
  % If indeed Ku equals 1/Kp, see Daubechies and Sweldens, 1998, Delft
  % book p 145
  % disp('Extra lifting steps, may use CDFI instead')
  Pac=Pa; Uac=Ua;
  clear Pa Ua
  Pa{1}=Pac;
  Ua{1}=Uac;
  Pa{2}=-1;
  Pa{3}=Ku;
  Ua{2}=Kp-1;
  Ua{3}=Kp-Kp^2;
  % Verify that this whole thing works here:
  difer([Kp 0 ; 0 Ku]-[1 Ua{3} ; 0 1]*[1 0 ; -Pa{3} 1]*...
	[1 Ua{2} ; 0 1]*[1 0 ; -Pa{2} 1]);
  Ku=1;
  Kp=1;
end

% Only for interior points
% SPLIT -----------------------------------------
d=x(even(x));
a=x(~even(x));

% Number of lifting steps------------------------
M=1;
if iscell(Pa); M=length(Pa); end

% FORWARD TRANSFORM
% Loop over lifting steps------------------------
for index=1:M
  if iscell(Pa)
    P=Pa{index};
    U=Ua{index};
  else
    P=Pa;
    U=Ua;
  end
  Pl=length(P);
  Ul=length(U);

  % PREDICT ---------------------------------------
  for l=ceil(Pl/2):ceil(length(x)/2)-floor(Pl/2)-(Pl==1)*mod(length(x),2)
    Lp=l+[(1-ceil(Pl/2)):1:floor(Pl/2)];
    if intel==1; d(l)=d(l)-floor(P(:)'*a(Lp)+1/2); end
    if intel==0; d(l)=d(l)-P(:)'*a(Lp); end
  end

  % UPDATE ----------------------------------------
  for l=1+floor(Ul/2):floor(length(x)/2)-ceil(Ul/2)+1
    Lu=l-[floor(Ul/2):-1:(1-ceil(Ul/2))];
    if intel==1; a(l)=a(l)+floor(U(:)'*d(Lu)+1/2); end
    if intel==0; a(l)=a(l)+U(:)'*d(Lu); end
  end
end

% SCALE -----------------------------------------
d=d*Ku;
a=a*Kp;

% Try the inverse transform by running it backwards
% UNDO SCALE ------------------------------------
y0=a/Kp;
y1=d/Ku;

if nargout>=2
  % INVERSE TRANSFORM
  % Reverse loop over lifting steps-------------------
  for index=M:-1:1
    if iscell(Pa)
      P=Pa{index};
      U=Ua{index};
    else
      P=Pa;
      U=Ua;
    end
    Pl=length(P);
    Ul=length(U);
    
    % UNDO UPDATE -----------------------------------
%    floor(length(x)/2)
    for l=1+floor(Ul/2):floor(length(x)/2)-ceil(Ul/2)+1
      Lu=l-[floor(Ul/2):-1:(1-ceil(Ul/2))];
      if intel==1; y0(l)=y0(l)-floor(U(:)'*y1(Lu)+1/2); end
      if intel==0; y0(l)=y0(l)-U(:)'*y1(Lu); end
    end
    
    % UNDO PREDICT ----------------------------------
%    ceil(length(x)/2)
%  ceil(Pl/2):ceil(length(x)/2)-floor(Pl/2)-(Pl==1)*mod(length(x),2)
    for l=ceil(Pl/2):ceil(length(x)/2)-floor(Pl/2)-(Pl==1)*mod(length(x),2)
      Lp=l+[(1-ceil(Pl/2)):floor(Pl/2)];
      if intel==1; y1(l)=y1(l)+floor(P(:)'*y0(Lp)+1/2); end
      if intel==0; y1(l)=y1(l)+P(:)'*y0(Lp); end
    end  
  end
  % MERGE -----------------------------------------
  xout(~even(x),:)=y0;
  xout(even(x),:)=y1;
  
  err=sqrt(mean((xout(:)-x(:)).^2));
  
  if err>1e-10
    warning(sprintf('Reconstruction failed with %8.2e',err))
  else
%    disp(sprintf('WTLIFT Reconstruction error %8.4e',err))
  end
%  disp(sprintf('WTLIFT (Analysis & Synthesis) took %8.4f s',...
%	       etime(clock,t0)))
end

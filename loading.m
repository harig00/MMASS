function [Hs,Gb,Zbav,C2bav,k,Zb,C2b,f,fK]=loading(R,T,Te,AT,kx,ky,F,method)
% [Hs,Gb,Zbav,C2bav,k,Zb,C2b,f,fK]=LOADING(R,T,Te,AT,kx,ky,F,method)
%
% For a model with an arbitrary number of layers
% (top layer is water or air and bottom layer is mantle halfspace),
% of density 'R' and depth 'T' and a loading scenario of
% load fractions 'F' and an elastic thickness 'Te', and for
% initial loads 'I', calculates a surface topography 'Hs'
% and a Bouguer gravity synthetic 'Gb'.
%
% INPUT
%
% R          density model for N layers and N-1 interfaces [kg/m3]
% T          depth to N-2 interfaces (first interface is "surface") [km]
% Te         elastic thickness [km]
% AT         applied topography at interface (dimension mXnX(N-1)) [km]
% kx, ky     wavenumber axes in rad/km (they better come out of KNUM2)
% F          load fractions - specify if AT is empty
% method     1 Using average loading ratio, not actual topography
%            2 Using frequency-dependent loading ratio calculated from loads
%            3 Using the actual loads if specified
%
%
% OUTPUT
%
% Hs         final topography in km
% Gb         Bouguer anomaly in mgal
% Zbav       Bouguer admittance in mgal/m (radially averaged)
% C2bav      Bouguer coherence square function (radially averaged)
%
% See MCKENZIE, LOADING1, KNUM2
%
% Last modified by fjsimons-at-alum.mit.edu, October 22nd, 2003

defval('R',[0 2600 2900 3300])
defval('T',[     0   15   35])
defval('Te',0)
defval('F',[     0.5 0.5   0])

defval('E',1.78e11); % Young's modulus
defval('v',0.25);    % Poisson's ratio
defval('method',3)   % See down

disp(sprintf('E= %5.3g; v= %5.3f; Te= %2.2i',E,v,Te))
disp(sprintf([ 'R= ' repmat('%6i',1,length(R))],R))
disp(sprintf([ 'T= ' repmat('%4i',1,length(T))],T))
g=fralmanac('GravAcc');
G=fralmanac('GravCst');

% Calculate wavenumbers ready for FFTSHIFT
% They have to come out of KNUM2 not FFTAXIS
[KX,KY]=meshgrid(kx(:),ky(:));
K=sqrt(KX.^2+KY.^2);

% Standard units of thickness in m
T=T*1000;
Te=Te*1000;
AT=AT*1000;
K=K/1000;

% Flexural Rigidity [Pa*m^3 or N*m]
D=(E*Te.^3)/(12*(1-v^2));

% Density contrasts: Number of interfaces
DR=diff(R);

% Some constants
K4=K.^4;
TpG=2*pi*G;
DK4g=D*K4/g;

% Initialize matrices
[ATK,ekz,Drekz,Z,Y,FK]=deal(repmat(NaN,size(AT)));
% Fourier transforms of demeaned fields
for index=1:length(T)
  if ~isempty(AT)
    % You specified the subsurface loading directly by AT
    % Long wavelengths in center, like the wavenumber matrix
    % Take all loading topographies and transform them to the
    % frequency domain
    ATK(:,:,index)=fftshift(fft2(AT(:,:,index)...
				 -mean(mean(AT(:,:,index)))));
  end
  % Upward continuation operators
  ekz(:,:,index)=exp(-K.*T(index));
  DRekz(:,:,index)=ekz(:,:,index)*DR(index);
end
% Calculate theoretical admittances of interface loading
Z(:,:,1)=-TpG*DR(1)./(DK4g+sum(DR(2:end))).*sum(DRekz(:,:,2:end),3);
for index=2:length(DR)
  inds=2:length(DR);
  inds=inds(inds~=index);
  jnds=1:length(DR);
  jnds=jnds(jnds~=index);
  Z(:,:,index)=TpG*(sum(DRekz(:,:,inds),3)-ekz(:,:,index).*...
		    (DK4g+sum(DR(jnds))));
end

% Convert to mgal/m
Z=Z/1e-5;

% Calculate Y's: turns LOADS (Pa) into TOPOGRAPHY (m); in 
Y(:,:,1)=(DK4g+sum(DR(2:end)))./(DK4g+sum(DR))/g/DR(1);
% Note that all subsequent Y's are identical
for index=2:length(DR)
  Y(:,:,index)=-1./(DK4g+sum(DR))/g;
end

% If AT is specified will calculate the loading (ratio) from the input
% data and not from the specified F
if ~isempty(AT)
  % Only the second interface is considered in this calculation
  fK=DR(2)*abs(ATK(:,:,2))/DR(1)./abs(ATK(:,:,1));
  fstd=std(real(fK(:)));
  f=mean(real(fK(:)));
  disp(sprintf('Calculated mean f= %5.3f Std f= %5.3f',f,fstd))
  F2=f/(1+f);
  F=[1-F2 F2 zeros([length(DR)-2 1])];
  for index=1:length(DR)
    FK(:,:,index)=DR(index)*g*abs(ATK(:,:,index));
  end  
  TL=repmat(sum(FK,3),[1 1 size(FK,3)]);
  FK=FK./TL;
end

% Figure out which program calls you and decide on a method
[p,n,e,v]=star69;
if strcmp(n,'loading1')
  method=1;
  elseif strcmp(n,'loading2')
  method=3;
end

% Now calculate topography
switch method  
 case 1
  % Without reference to the actual topography; hence
  % purely theoretical result
  HK=shiftdim(repmat(F(:),[1 size(DK4g)]),1).*Y;
 case 2
  % Uses frequency-dependent loading ratio calculated from
  % the absolute value of the data but doesn't use the actual
  % data itself anymore
  HK=FK.*Y;
 case 3
  % Works with the actual data input
  % ATK is in m
  % Y is in 1/rho/g
  HK=ATK.*Y.*shiftdim(repmat(DR(:),[1 size(DK4g)]),1)*g;
end

% Message
msg=[{ 'Using average loading ratio.'} ...
      { 'Using frequency-dependent loading ratio.'} ...
      { 'Using actual loads'}];
disp(msg{method})

% Calculate Bouguer gravity anomaly
GBK=HK.*Z;

% Calculate unaveraged "spectra"
% So because it's exact can just add it. For real
% data need the smoothing.
SHG=sum(GBK.*conj(HK),3);
HH=sum(HK.*conj(HK),3);
GG=sum(GBK.*conj(GBK),3);

% Calculate Bouguer admittance
Zb=SHG./HH;
% Calculate Bouguer coherence
C2b=abs(SHG).^2./HH./GG;

% Calculate radially averaged functions
[k,Zbav]=isav(Zb,kx,ky); Zbav=Zbav(:);
[k,C2bav]=isav(C2b,kx,ky); C2bav=C2bav(:);

% Calculate topography and gravity in the spatial domain
if ~isempty(AT)
      [Hs,Gb]=deal(repmat(0,size(DK4g)));
      % The way I do this now with KNUM2, IFFT2 should return reals
      for index=1:length(DR)
	Hsc=ifft2(ifftshift(HK(:,:,index)));
	Gbc=ifft2(ifftshift(GBK(:,:,index)));
	if mean(abs(imag(Hsc(:))))<1e-14
	  Hs=Hs+real(Hsc);
	else 
	  warning(sprintf('Transformed topography is not real by %s',...
			  mean(abs(imag(Hsc(:))))))
	end
	if mean(abs(imag(Gbc(:))))<1e-14
	  Gb=Gb+real(Gbc);
	else 
	  warning(sprintf('Transformed gravity is not real by %s',...
			  mean(abs(imag(Gbc(:))))))
	end
      end
      Hs=Hs/1000; % In km
else
  [Hs,Gb]=deal(NaN);
end



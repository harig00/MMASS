function [l,Zb,G2b,Zf,G2f]=mckenzie(R,T,F,Te)
% [l,Zb,G2b,Zf,G2f]=MCKENZIE(R,T,F,Te)
%
% This program works without specifying the actual topography!
%
% For a model with an arbitrary number of layers
% (top layer is water or air and bottom layer is mantle halfspace),
% calculates free-air and Bouguer admittances and coherences
% as a function of wavelength, for a given flexural rigidity.
% Note that for a two-layer model, the depth of the compensating
% interface does not matter at all for the coherence.
%
% INPUT:
%
% R       density model for N layers and N-1 interfaces in kg/m3
% T       depth to N-2 interfaces (first interface is "surface") in km
% F       fraction of total load at each of N-1 interface
% Te      elastic thickness in km, cannot be a vector
%
% OUTPUT:
%
% l       wavelength in m
% Zb      Bouguer admittance in mgal/m
% G2b     Bouguer coherence-square
% Zf      Free-air admittance in mgal/m
% G2f     Free-air coherence-square
%
% EXAMPLE: Forsyth's Fig. 1
%
% Te=[40 30 20 10 5 0];
% [l,Zb,G2b]=mckenzie([0 2670 2670+670],[0 35],[1 0],Te(4)); % Top loading
% [l,Zb,G2b]=mckenzie([0 2670 2670+670],[0 35],[0 1],Te(4)); % Bottom loading
%
% See FORSYTH, MCKENZIE1, MCKENZIE2, LOADING
%
% Last modified by fjsimons-at-alum.mit.edu, October 21st, 2003

defval('R',[1030 2600 2900 3300])
defval('T',[      0     15   35])
defval('F',[    0.6    0.4    0])
defval('Te',17)
defval('E',1.4e11);
v=0.25;
g=9.81;
disp(sprintf('E= %5.3g; v= %5.3f',E,v))
G=fralmanac('GravCst');
l=linspace(1,3000,2000)*1000;
k=2*pi./l;
k4=k.^4;

% Standard units
T=T*1000;
Te=Te*1000;

% Flexural Rigidity [Pa*m^3 or N*m]
D=(E*Te.^3)/(12*(1-v^2));

TpG=2*pi*G;
Dk4g=D*k4/g;

DR=diff(R);

[K,ZS]=meshgrid(k,T);
ekz=exp(-K.*ZS);
DRekz=ekz.*repmat(DR(:),1,size(K,2));

% Calculate Z's
Z(1,:)=-TpG*DR(1)./(Dk4g+sum(DR(2:end))).*sum(DRekz(2:end,:),1);
for index=2:size(DR(:),1)
  inds=2:length(DR);
  inds=inds(inds~=index);
  jnds=1:length(DR);
  jnds=jnds(jnds~=index);
  Z(index,:)=TpG*(sum(DRekz(inds,:),1)-ekz(index,:).*(Dk4g+sum(DR(jnds))));
end

% Convert to mgal/m
Z=Z/1e-5;

% Calculate Y's
Y(1,:)=(Dk4g+sum(DR(2:end)))./(Dk4g+sum(DR))/g/DR(1);
for index=2:size(DR(:),1)
  Y(index,:)=-1./(Dk4g+sum(DR))/g;
end

% Calculate Bouguer admittance
Zb=((F(:)'.^2)*(Y.^2.*Z))./((F(:)'.^2)*(Y.^2));

% Calculate Bouguer coherence-square
G2b=(([F(:)'].^2)*(Y.^2.*Z)).*(([F(:)'].^2)*(Y.^2.*Z))./...
   (([F(:)'].^2)*(Y.^2))./(([F(:)'].^2)*(Y.^2.*Z.^2));

% Calculate Free-air admittance in mgal/m
Z=Z+repmat(TpG*DR(1),size(Z))/1e-5;
Zf=((F(:)'.^2)*(Y.^2.*Z))./((F(:)'.^2)*(Y.^2));

% Calculate Free-air coherence-square
G2f=(([F(:)'].^2)*(Y.^2.*Z)).*(([F(:)'].^2)*(Y.^2.*Z))./...
   (([F(:)'].^2)*(Y.^2))./(([F(:)'].^2)*(Y.^2.*Z.^2));

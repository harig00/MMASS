function [Sstrain,Tstrain]=modestrain(rs,re,rfs,rad,nn,el,ww,U,V,P,dUdr,dVdr,W,dWdr)
% [Sstrain,Tstrain]=MODESTRAIN(rs,re,rfs,rad,nn,el,ww,U,V,P,dUdr,dVdr,W,dWdr)
%
% Gets coefficients of the strain as pertains to modal excitation.
% The input comes out of GETSPHEROIDAL and GETTOROIDAL.
%
% INPUT:
%
% rs          Radius at which strain functions are to be returned (e.g. source) [m]
% re          Radius at which P is to be returned (e.g. surface) [m]
% rfs         Radius at which U is to be returned (e.g. solid surface) [m]
% rad         Radius at which the input radial eigenfunctions are evaluated [m]
% nn          Array of integer branch numbers
% el          Array of spherical harmonic degrees
% ww          Array of angular frequencies [rad/s]
% U, dUdr     Radial-displacement eigenfunctions and derivatives [kg^{-1/2}, kg^{-1/2}/m]
% V, dVdr     Tangential-displacement eigenfunctions, derivatives [kg^{-1/2}, kg^{-1/2}/m]
% P           Potential-perturbation eigenfunctions [units, units/m]
% W, dWdr     Radial radial eigenfunctions and their derivatives [kg^{-1/2}, kg^{-1/2}/m]
%
% OUTPUT
%
% Sstrain     A structure array with nn,el,ww and these eigenfunction combinations:
%
%   a         The functions multiplying the Mrr monopole source at rs [units/m]
%   b         The functions multiplying the quadrupole sources at rs [units/m]
%   c         The functions multiplying the Mtt+Mpp monopole source at rs [units/m]
%   d         The functions multiplying the dipole sources at rs [units/m]
%   Pre       Potential-perturbation radial eigenfunctions at re [units]
%   Urfs      Radial-displacement radial eigenfunctions at rfs [units]
%   rs        The actual radius rs used, as returned from FINDRAD [m]
%   re        The actual radius rs used, as returned from FINDRAD [m]
%   rfs       The actual radius rs used, as returned from FINDRAD [m]
%
% Tstrain     A structure array with nn,el,ww and these eigenfunction combinations:
%  
%   bla
%   di
%   bla
%
% Last modified by fjsimons-at-alum.mit.edu, 05/31/2010

% Check the radii are sorted so that the "end" is the Earth surface
% (ocean or not is immaterial for the potential that we return here)
difer(rad-sort(rad),[],[],NaN)

defval('re',rad(end))
% This default presupposed it's PREM with an ocean.
defval('rfs',6368000)

% Evaluate the strain functions at the requested radius, e.g. the source
rsint=findrad(rad,rs); rs=rad(rsint);

% Evaluate the potential eigenfunctions at the requested radius
reint=findrad(rad,re,2); re=rad(reint);
% Evaluate the displacement eigenfunctions - remember FINDRAD returns the upperside
rfsint=findrad(rad,rfs,2);
if ~difer(rad(rfsint)-rad(rfsint-1))
  rfsint=rfsint-1;
  disp('Taking the underside rather than the upperside of this discontinuity')
end
rfs=rad(rfsint);

% Restrict immediately to the requested radii
Urfs=U(rfsint,:);
Pre=P(reint,:);
U=U(rsint,:);
V=V(rsint,:);
dUdr=dUdr(rsint,:);
dVdr=dVdr(rsint,:);

% Produce the relevant combinations
% These guys all have the same units per m
a=dUdr;
b=1/2*V/rs;
c=1/2*(2*U-repmat(el(:)'.*(el(:)'+1),size(V,1),1).*V)/rs;
d=dVdr+(U-V)/rs;

% Produce output in a structure array
Sstrain=struct('nn',nn,'el',el,'ww',ww,...
	      'a',a,'b',b,'c',c,'d',d,...
	       'Pre',Pre,'Urfs',Urfs,...
	       'rs',rs,'re',re,'rfs',rfs);

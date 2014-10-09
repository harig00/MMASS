function Sstrain=smodestrain(rs,re,rfs,rad,nn,el,ww,U,V,P,dUdr,dVdr)
% Sstrain=SMODESTRAIN(rs,re,rfs,rad,nn,el,ww,U,V,P,dUdr,dVdr)
%
% Gets coefficients of the strain from spheroidal-mode excitation. 
%
% INPUT:
%
% rs          Radius at which strain functions are to be returned [m],
%             typically at the source, where the excitation is computed
% re          Radius at which P is to be returned [m],
%             somewhere else, e.g. at the ocean surface or at satellite height
% rfs         Radius at which U is to be returned [m],
%             typically somewhere else, e.g. at the fluid-solid interface.
%             And then follow the modal eigenfunctions from GETSPHEROIDAL:
% rad         Radius at which the input radial eigenfunctions are evaluated [m]
% nn          Array of integer branch numbers
% el          Array of spherical harmonic degrees
% ww          Array of angular frequencies [rad/s]
% U, dUdr     Radial-displacement eigenfunctions and derivatives 
%             Although "displacement", units are [kg^{-1/2}, kg^{-1/2}/m]
% V, dVdr     Tangential-displacement eigenfunctions, derivatives
%             Although "displacement", units are [kg^{-1/2}, kg^{-1/2}/m]
% P           Potential-perturbation radial eigenfunctions
%             Although "potential", units are [m/s^2/sqrt(kg) 1/s^2/sqrt(kg)]
%
% OUTPUT
%
% Sstrain     A structure array with nn,el,ww and these eigenfunction combinations:
%
%   a         The functions multiplying the Mrr monopole source at rs
%   b         The functions multiplying the quadrupole sources at rs 
%   c         The functions multiplying the Mtt+Mpp monopole source at rs
%   d         The functions multiplying the dipole sources at rs
%             Although "strain", units are [kg^{-1/2}/m]
%   Pre       Potential-perturbation radial eigenfunctions at re
%   Prfs      Potential-perturbation radial eigenfunctions at rfs
%   Urfs      Radial-displacement radial eigenfunction U at rfs
%   Vrfs      Radial-displacement radial eigenfunction V at rfs
%   rs        The actual radius rs used, as returned from FINDRAD [m]
%   re        The actual radius re used, as returned from FINDRAD [m]
%   rfs       The actual radius rfs used, as returned from FINDRAD [m]
%
% SEE ALSO: 
%
% GETSPHEROIDAL, TMODESTRAIN
%
% Last modified by fjsimons-at-alum.mit.edu, 07/21/2010

defval('xver',0)

% Check the radii are sorted so that the "end" is the Earth surface
% (ocean or not is immaterial for the potential that we return here)
difer(rad-sort(rad),[],[],NaN)

defval('re',rad(end))
% This default presupposes it's PREM with an ocean.
defval('rfs',6368000)

% Evaluate the strain functions at the requested radius, e.g. the source
rsint=findrad(rad,rs,1*[xver~=0]+3*[xver==0]); rs=rad(rsint);

% Evaluate the potential eigenfunctions at the other requested radius
reint=findrad(rad,re,2*[xver~=0]+3*[xver==0]); re=rad(reint);

% Evaluate the displacement eigenfunctions - remember FINDRAD returns the upperside
rfsint=findrad(rad,rfs,2*[xver~=0]+3*[xver==0]);
if ~difer(rad(rfsint)-rad(rfsint-1))
  rfsint=rfsint-1;
  disp('Taking the underside rather than the upperside of this discontinuity')
end
rfs=rad(rfsint);

% Restrict immediately to the requested radii
% This is for the point at which we want to know the displacement result
Urfs=U(rfsint,:);
Vrfs=V(rfsint,:);
Prfs=P(rfsint,:);

% So here later should build in upward continuation by satellite
Pre=P(reint,:);

% This is for the excitation, after which I forget it
 U  = U(rsint,:);
dUdr=dUdr(rsint,:);
 V  = V(rsint,:);
dVdr=dVdr(rsint,:);

% Produce the relevant combinations by DT p. 371
% These guys all have the same units of [1/m/sqrt(kg)]
a=dUdr;
b=1/2*V/rs;
c=1/2*(2*U-repmat(el(:)'.*(el(:)'+1),size(V,1),1).*V)/rs;
d=dVdr+(U-V)/rs;

% Produce output in a structure array
Sstrain=struct('nn',nn,'el',el,'ww',ww,...
	      'a',a,'b',b,'c',c,'d',d,...
	       'Pre',Pre,...
	       'Prfs',Prfs,'Urfs',Urfs,'Vrfs',Vrfs,...
	       'rs',rs,'re',re,'rfs',rfs);

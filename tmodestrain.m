function Tstrain=tmodestrain(rs,rfs,rad,nn,el,ww,W,dWdr)
% Tstrain=TMODESTRAIN(rs,rfs,rad,nn,el,ww,W,dWdr) 
%
% Gets coefficients of the strain from toroidal-mode excitation.
%
% INPUT:
%
% rs          Radius at which strain functions are to be returned [m],
%             typically at the source, where the excitation is computed
% rfs         Radius at which W is to be returned [m],
%             typically somewhere else, e.g. at the fluid-solid interface.
%             And then follow the modal eigenfunctions from GETTOROIDAL:
% rad         Radius at which the input radial eigenfunctions are evaluated [m]
% nn          Array of integer branch numbers
% el          Array of spherical harmonic degrees
% ww          Array of angular frequencies [rad/s]
% W, dWdr     Toroidal radial eigenfunctions and derivatives
%             Although "displacement", units are [kg^{-1/2}, kg^{-1/2}/m]
%
% OUTPUT:
%
% Tstrain     A structure array with nn,el,ww and these eigenfunction combinations:
%
%   e         The functions multiplying the quadrupole sources at rs
%   f         The functions multiplying the dipole sources at rs
%             Although "strain", units are [kg^{-1/2}/m]
%   Wrfs      Radial-displacement radial eigenfunction U at rfs
%   rs        The actual radius rs used, as returned from FINDRAD [m]
%   rfs       The actual radius rfs used, as returned from FINDRAD [m]
%
% SEE ALSO: 
%
% GETTOROIDAL, SMODESTRAIN
%
% Last modified by efwelch@princeton.edu 07/16/2010
% Last modified by fjsimons-at-alum.mit.edu, 07/16/2010

% Check the radii are sorted so that the "end" is the Earth surface
% (ocean or not is immaterial for the potential that we return here)
difer(rad-sort(rad),[],[],NaN)

% This default presupposes it's PREM with an ocean.
defval('rfs',6368000)

% Evaluate the strain functions at the requested radius, e.g. the source
rsint=findrad(rad,rs); rs=rad(rsint);

% Evaluate the displacement eigenfunctions - remember FINDRAD returns the upperside
rfsint=findrad(rad,rfs,2);
if ~difer(rad(rfsint)-rad(rfsint-1))
  rfsint=rfsint-1;
  disp('Taking the underside rather than the upperside of this discontinuity')
end
rfs=rad(rfsint);

% Restrict immediately to the requested radii
% This is for the point at which we want to know the displacement result
Wrfs=W(rfsint,:);
% This is for the excitation, after which I forget it
 W  = W(rsint,:);
dWdr=dWdr(rsint,:);

% Produce the relevant combinations by DT p. 371
% These guys all have the same units per m
e=(dWdr-W./rs);
f=W/rs;

% Produce output in a structure array
Tstrain=struct('nn',nn,'el',el,'ww',ww,...
	      'e',e,'f',f,...
	       'Wrfs',Wrfs,...
	       'rs',rs,'rfs',rfs);

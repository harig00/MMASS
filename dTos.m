function varargout=dTos(k,th,params,phi,xi,pxm)
% [dTinvdD,dTinvdf2]=dTos(k,th,params,phi,xi,pxm)
%
% Calculates first derivatives of isotropic Tinv in Olhede & Simons
% under the UNCORRELATED equilibrium-loading model
%
% k          Wavenumber(s) at which this is to be evaluated [1/m]
% th         The parameter vector with TWO elements (rest ignored)
%            th(1)=D    Isotropic flexural rigidity 
%            th(2)=f2   The sub-surface to surface initial loading ratio 
% params     A structure with AT LEAST these constants that are known:
%            DEL   surface and subsurface density contrast [kg/m^3]
%            g     gravitational acceleration [m/s^2]
% phi        Optionally precalculated phi, see PHIOS
% xi         Optionally precalculated xi, see XIOS
% pxm        Optionally precalculated (phi*xi-1)
%
% OUTPUT:
%
% dTinvdD    The first derivative of Tinv with respect to D
% dTinvdf2   The first derivative of Tinv with respect to f2
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2012

% Extract the parameters from the input
D=th(1);
f2=th(2);
DEL=params.DEL;
g=params.g;

if nargin<6
  % First the auxiliary quantities
  phi=phios(k,D,DEL,g);
  xi = xios(k,D,DEL,g);
  % Note that this has a zero at zero wavenumber
  pxm=(phi.*xi-1);
end

% First dTinvdD, which is a symmetric matrix
warning off MATLAB:divideByZero
fax=-2*k(:).^4/g/f2*dpos(DEL,-2,1).*xi./pxm.^3;
warning on MATLAB:divideByZero
dTinvdD(:,1)=fax.*(1+f2*dpos(DEL,2,-2)*phi.^2+f2*dpos(DEL,1,-1)*phi.*xi+... 
		   dpos(DEL,-1,1).*xi.^2);  
dTinvdD(:,2)=fax.*(f2*dpos(DEL,1,-1)*phi+1/2*f2*xi+3/2*dpos(DEL,-1,1)*xi+...
		   1/2*f2*phi.*xi.^2-1/2*dpos(DEL,-1,1)*phi.*xi.^2+...
		   dpos(DEL,-2,2)*xi.^3);
dTinvdD(:,3)=fax.*(f2+f2*dpos(DEL,-1,1)*xi.^2+2*dpos(DEL,-2,2)*xi.^2-...
		   dpos(DEL,-2,2)*phi.*xi.^3+dpos(DEL,-3,3)*xi.^4);

if nargout>=2
  % Then dTinvdf2, which is also a symmetric matrix
  warning off MATLAB:divideByZero
  fax=-dpos(DEL,-2,2)*xi.^2./pxm.^2/f2^2;
  warning on MATLAB:divideByZero
  dTinvdf2(:,1)=fax;
  dTinvdf2(:,2)=fax.*dpos(DEL,-1,1).*xi;
  dTinvdf2(:,3)=fax.*dpos(DEL,-2,2).*xi.^2;
else
  dTinvdf2=NaN;
end

% Output
varns={dTinvdD,dTinvdf2};
varargout=varns(1:nargout);

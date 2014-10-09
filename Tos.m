function varargout=Tos(k,th,params,phi,xi,pxm,xver)
% [invT,detT,L,T]=Tos(k,th,params,phi,xi,pxm,xver)
%
% Calculates functions of isotropic T in Olhede & Simons for the
% UNCORRELATED loading scenario with the primary spectra S11 being the
% equilibrium-loading ones. 
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with TWO elements (rest ignored)
%          D   Isotropic flexural rigidity [Nm]
%          f2  The sub-surface to surface initial loading ratio
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% phi      Optionally precalculated phi, see PHIOS
% xi       Optionally precalculated xi, see XIOS
% pxm      Optionally precalculated (phi*xi-1)
%
% OUTPUT:
%
% invT     A 3-column vector with all the wavenumbers unwrapped,
%          invT={invT[1,1](k) invT[1,2](k) invT[2,2](k)}
% detT     The determinant of T, a column vector over the wavenumbers
% L        The Cholesky factorization of T, as the lower-left matrix 
%          L={L[1,1](k) L[2,1](k) L[2,2](k)}, where L[1,2]=0
% T        The actual matrix T, in the same format as invT
%
% EXAMPLE:
%
% [~,~,th0,p,k]=simulos([],[],[],1);
% [invTo,detTo,Lo,To]=Tos(k,th0,p,[],[],[],1);
%
% Last modified by fjsimons-at-alum.mit.edu, 02/10/2011

defval('xver',0)

% Extract the parameters from the input
D=th(1);
f2=th(2);
DEL=params.DEL;
g=params.g;

defval('phi',phios(k,D,DEL,g));
defval('xi',xios(k,D,DEL,g));
% Note that this has a zero at zero wavenumber
defval('pxm',(phi.*xi-1));

% Forcefully set f2 to a positive number even if it means a throw back
f2=abs(f2);

% The inverse of T; ignore warnings as Inf turns to NaN in HFORMOS 
warning off MATLAB:divideByZero
fax=dpos(DEL,-2,2)*xi.^2/f2./pxm.^2;
warning on MATLAB:divideByZero
invT=[fax.*(                   1+f2*dpos(DEL,2,-2)*phi.^2) ...
      fax.*(dpos(DEL,-1,1)*xi   +f2*dpos(DEL,1,-1)*phi   ) ...
      fax.*(dpos(DEL,-2,2)*xi.^2+f2)];

if nargout>=2 || xver==1
  % Compute the determinant of T; this will be zero at k=0
  detT=f2*dpos(DEL,4,-4)*xi.^(-4).*pxm.^2;
else
  detT=NaN;
end

if nargout>=3 || xver==1
  % Compute the Cholesky factorization of T
  fax=dpos(DEL,0,-1)*xi.^(-1)./...
      sqrt(dpos(DEL,0,2)*xi.^2+f2*dpos(DEL,2,0));
  L=[fax.*( dpos(DEL,0,2)*xi.^2+f2*dpos(DEL,2,0)) ...
     fax.*(-dpos(DEL,1,-1)*[dpos(DEL,0,2)*xi+f2*dpos(DEL,2,0).*phi]) ...
     fax.*(sqrt(f2)*dpos(DEL,2,0)*pxm)];
else
  L=NaN;
end

if nargout>=4 || xver==1
  % Compute T itself, which is required when producing blurred things
  fax=xi.^(-2);
  T=[fax.*(             xi.^2+f2*dpos(DEL,2,-2)        ) ...
     fax.*(-dpos(DEL,1,-1)*xi-f2*dpos(DEL,3,-3)*phi    ) ...
     fax.*( dpos(DEL,2,-2)   +f2*dpos(DEL,4,-4)*phi.^2)];
else
  T=NaN;
end

% And now for the output
varns={invT,detT,L,T};
varargout=varns(1:nargout);

% Verification mode
if xver==1
  disp('Tos being verified')
  % Explicit verification of the determinant
  detcheck(detT ,T ,8)

  % Explicit verification of the inverse by the Cayley-Hamilton theorem
  invcheck(invT ,detT ,T ,3,1)

  % Explicit verification of the inverse by checking the identity
  invcheck(invT ,detT ,T ,3,2)

  % Check the Cholesky by multiplication
  cholcheck(L,T,5,1)

  % Check the Cholesky by factorization at a random wave number
  cholcheck(L,T,4,2)
end

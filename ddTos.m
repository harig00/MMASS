function varargout=ddTos(k,th,params,phi,xi,pxm,xver)
% [trWfax,dTinvdDdf2]=ddTos(k,th,params,phi,xi,pxm,xver)
%
% Calculates functions of the second derivatives of isotropic Tinv in
% Olhede & Simons (2013). 
%
% k           Wavenumber(s) at which this is to be evaluated [1/m]
% th          The parameter vector with TWO elements (rest ignored)
%             D   Isotropic flexural rigidity [Nm]
%             f2  The sub-surface to surface initial loading ratio
% params      A structure with AT LEAST these constants that are known:
%             DEL   surface and subsurface density contrast [kg/m^3]
%             g     gravitational acceleration [m/s^2]
% phi         Optionally precalculated phi, see PHIOS
% xi          Optionally precalculated xi, see XIOS
% pxm         Optionally precalculated (phi*xi-1)
%
% OUTPUT:
%
% trWfax      The sum of the eigenvalues of dTinvdDdf2
% dTinvdDdf2  The second (cross) derivative of Tinv to D and f2
%
% Last modified by fjsimons-at-alum.mit.edu, 12/20/2010

defval('xver',0)

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

% Produce the sum of the eigenvalues
warning off MATLAB:divideByZero
fax=2*k(:).^4/g/f2*dpos(DEL,-1,-1).*xi.^(-1)./pxm;
warning on MATLAB:divideByZero
trWfax=fax.*(dpos(DEL,0,1)*xi.^2-dpos(DEL,1,0)*pxm);

if nargout>=2 || xver==1
  % First dTinvdDdf2, which is a symmetric matrix
  warning off MATLAB:divideByZero
  fax=2*k(:).^4/g/f2^2*dpos(DEL,-2,1).*xi./pxm.^3;
  warning on MATLAB:divideByZero
  dTinvdDdf2(:,1)=fax.*(1+dpos(DEL,-1,1)*xi.^2);  
  dTinvdDdf2(:,2)=fax.*(dpos(DEL,-1,1)*...
			[dpos(DEL,-1,1)*xi.^3-1/2*phi.*xi.^2+3/2*xi]);
  dTinvdDdf2(:,3)=fax.*(dpos(DEL,-2,2)*...
			[dpos(DEL,-1,1)*xi.^4-phi.*xi.^3+2*xi.^2]);
else
  dTinvdDdf2=NaN;
end

if xver==1
  [invT,detT,L]=Tos(k,th,params,phi,xi,pxm);
  randi=ceil(rand*length(dTinvdDdf2));
  M=dTinvdDdf2;
  Mrand=[M(randi,1) M(randi,2) ; M(randi,2) M(randi,3)];
  Lrand=[L(randi,1)    0       ; L(randi,2) L(randi,3)];
  % Check RELATIVE precision
  difer((trWfax(randi)-trace(Lrand'*Mrand*Lrand))/trWfax(randi))
end

% And now for the output
varns={trWfax,dTinvdDdf2};
varargout=varns(1:nargout);

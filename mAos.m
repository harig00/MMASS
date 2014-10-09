function varargout=mAos(k,th,params,phi,xi,pxm,xver)
% [m,A]=mAos(k,th,params,phi,xi,pxm,xver)
%
% Calculates some auxiliary numbers and matrices for Olhede & Simons
% under the UNCORRELATED loading model
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with elements:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% phi      Optionally precalculated phi, see PHIOS
% xi       Optionally precalculated xi, see XIOS
% pxm      Optionally precalculated (phi*xi-1)
%
% OUTPUT:
%
% m        The "means" parameters which enter in the score, as a cell
% A        The "A" matrices which enter in the score, as a cell
%
% EXAMPLE:
% 
% [~,~,th0,p,k]=simulos([],[],[],1);
% [m,A]=mAos(k,th0,p,[],[],[],1);
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2012

defval('xver',0)

% Extract the parameters from the input
D=th(1);
f2=th(2);
s2=th(3);
nu=th(4);
rho=th(5);
DEL=params.DEL;
g=params.g;

% First the auxiliary quantities
defval('phi',phios(k,D,DEL,g));
defval('xi', xios(k,D,DEL,g));
% Note that this has a zero at zero wavenumber
defval('pxm',(phi.*xi-1))

% A variable that is also needed
avark=4*nu/pi^2/rho^2+k(:).^2;

% First calculate the "means" that depend on the wavenumbers
warning off MATLAB:divideByZero
m{1}=k(:).^4/g/dpos(DEL,0,1).*(2-phi.*xi+dpos(DEL,-1,1).*xi.^2)./xi./pxm;
warning on MATLAB:divideByZero
m{2}=1/2/f2;
m{3}=1/s2;
m{4}=(nu+1)/nu+log(4*nu/pi^2/rho^2)-...
     4*(nu+1)/pi^2/rho^2./avark-log(avark);
m{5}=-2*nu/rho+8*nu/rho*(nu+1)/pi^2/rho^2./avark;

% Full output and extra verification etc
if nargout>1 || xver==1
  % Initialize
  A=cellnan(length(th),length(k(:)),3);
  % Then calculate the "matrices" that depend on the wavenumbers
  [A{1},A{2}]=dTos(k,th,params,phi,xi,pxm);
  if xver==0
    % Get the inverse of T which enters the below
    invT=Tos(k,th,params,phi,xi,pxm);
  else
    % Get the inverse of T and the Cholesky for later testing
    [invT,detT,L]=Tos(k,th,params,phi,xi,pxm);
  end
  A{3}=-m{3}*invT;
  A{4}=-repmat(m{4},1,3).*invT;
  A{5}=-repmat(m{5},1,3).*invT;

  % Verification mode
  if xver==1
    % Check various identities at random wavenumbers 
    tracecheck(L,A,m,9,1)
  end
else
  A=NaN;
end

% Output
varns={m,A};
varargout=varns(1:nargout);

function mcF=Fisherkos(k,th,params,xver)
% mcF=Fisherkos(k,th,params,xver)
%
% Calculates the entries in the Fisher matrix for Olhede & Simons (2013)
% for the UNCORRELATED initial-loading model prior to wavenumber averaging. 
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The five-parameter vector with elements:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
%
% OUTPUT:
%
% mcF      The 15-column Fisher matrix, in this order:
%          [1]  FDD    [2]  Ff2f2
%          [3]  Fs2s2  [4]  Fnunu [5] Frhorho
%          [6]  FDf2   [7]  FDs2  [8]  FDnu  [9] FDrho
%          [10] Ff2s2  [11] Ff2nu [12] Ff2rho
%          [13] Fs2nu  [14] Fs2rho
%          [15] Fnurho
%
% EXAMPLE:
% 
% [~,~,th0,p,k]=simulos([],[],[],1);
% mcF=Fisherkos(k,th0,p,1);
%
% Last modified by fjsimons-at-alum.mit.edu, 02/10/2011

defval('xver',0)

% Extract the parameters from the input
D=th(1);
f2=th(2);
DEL=params.DEL;
g=params.g;

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% First the auxiliary quantities
phi=phios(k,D,DEL,g);
xi = xios(k,D,DEL,g);
% Note that this has a zero at zero wavenumber
pxm=(phi.*xi-1);

% First compute the "means" parameters
m=mAos(k,th,params,phi,xi,pxm,xver);

% Forcefully set f2 to a positive number even if it means a throw back
f2=abs(f2);

% Then compute the various entries in the order of the paper
lk=length(k(:));
% Some of them depend on the wave vectors, some don't
mcF=cellnan([npp 1],...
	    [lk 1 1 repmat(lk,1,6) 1  repmat(lk,1,5)],...
	    repmat(1,1,npp));

% FDD
warning off MATLAB:divideByZero
fax=2*k(:).^8/g^2*dpos(DEL,-2,-2)./pxm.^2.*xi.^(-2);
warning on MATLAB:divideByZero
mcF{1}=fax.*(2*dpos(DEL,0,2)*xi.^4+2*dpos(DEL,2,0)*phi.^2.*xi.^2+...
	     4*dpos(DEL,2,0)-4*dpos(DEL,1,1)*phi.*xi.^3+6*dpos(DEL,1,1)*xi.^2-... 
	     4*dpos(DEL,2,0)*phi.*xi+f2*dpos(DEL,2,0)*xi.^2+...
	     dpos(DEL,0,2)*xi.^2/f2);
% Ff2f2
mcF{2}=1/f2^2;

% Fthsths
for j=3:np
  mcF{j}=2*m{j}.^2;
end

% FDf2 might write it in explicitly, compare to FISHERKROS
mcF{np+1}=ddTos(k,th,params,phi,xi,pxm);

% All further combinations
jcombo=combntns(1:np,2);
for j=2:length(jcombo)
  mcF{np+j}=2*m{jcombo(j,1)}.*m{jcombo(j,2)};
end

% Verification step
if xver==1
  % The "linear" tracecheck has already happed in MAOS
  % This should be twice the negative sum of the eigenvalues of
  % the product of this...
  [M,N]=dTos(k,th,params,phi,xi,pxm);
  % and this...
  [invT,detT,L]=Tos(k,th,params,phi,xi,pxm);
  % ... which we'll be again checking at a random wavenumber
  tracecheck(L,{M N},{m{1:2}},9,1)

  % ... and then these also for ueber-completeness
  for ind=3:np
    tracecheck(L,{-repmat(m{ind},[size(invT,1)/size(m{ind},1)],3).*invT},...
	       {m{ind}},9,1)
  end

  % And then we check the sum of the squares of the eigenvalues which
  % hasn't happened so far
  tracecheck(L,{M N},{mcF{1:2}},8,2)
end

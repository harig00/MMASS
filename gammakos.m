function gammak=gammakos(k,th,params,Hk)
% gammak=GAMMAKOS(k,th,params,Hk)
%
% Computes the first partial derivatives of the likelihood function for
% the UNCORRELATED isotropic Forsyth model in Olhede & Simons. 
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The five-parameter vector argument [scaled]:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% Hk       A complex matrix of Fourier-domain observations
%
% OUTPUT:
%
% gammak   A 5-column vector with the wavenumbers unwrapped, containing
%          the five partials of the likelihood function as columns.
%
% Last modified by fjsimons-at-alum.mit.edu, 07/08/2013

% Extract the needed parameters from the input
D=th(1);
DEL=params.DEL;
g=params.g;

% The number of parameters to solve for
np=5;

% First the auxiliary quantities
phi=phios(k,D,DEL,g);
xi = xios(k,D,DEL,g);
% Note that this has a zero at zero wavenumber
pxm=(phi.*xi-1);

% First get the special matrices etc.
[m,A]=mAos(k,th,params,phi,xi,pxm);

% Then get the power spectrum
S11=maternos(k,th);

% Now compute the score properly speaking
gammak=nan(length(k(:)),np);
for j=1:np
  gammak(:,j)=-2*m{j}-hformos(S11,A{j},Hk);
end


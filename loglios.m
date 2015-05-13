function [L,gam,momx]=loglios(th,params,Hk,k,scl)
% [L,gam,momx]=LOGLIOS(th,params,Hk,k,scl)
%
% Calculates the full negative logarithmic likelihood and its
% derivatives, i.e. minus LKOS and minus GAMMAKOS averaged over
% wavenumber space. This is the function that we need to MINIMIZE!
%
% INPUT:
%
% th       The five-parameter vector argument [scaled]:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%           kiso   wavenumber beyond which we are not considering the likelihood
% Hk       A [prod(params.NyNx)*2]-column vector of complex Fourier-domain observations
% k        the wavenumbers at which these are being evaluated [1/m]
% scl      The vector with any scalings applied to the parameter vector
% momx     Moments of the quadratic piece Xk over all relevant wavenumbers
%
% OUTPUT:
%
% L        The loglihood, averaged over all wavenumbers
% gam      The score, averaged over all wavenumbers
%
% SEE ALSO:
%
% FISHERKOS, which should be incorporated at a later stage
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2014

% Default scaling is none
defval('scl',ones(size(th)))

% Scale up the parameter vector for the proper likelihood and score
th=th.*scl;

% Here I build the protection that the flexural rigidity,
% subsurface-to-surface ration, and the three Matern parameters should be
% positive. I mirror them up! Thereby messing with the iteration path,
% but hey. It means we can use FMINUNC also.
th([1 2 3 4 5])=abs(th([1 2 3 4 5]));

% Filter, perhaps
[Lk,Xk]=Lkos(k,th,params,Hk);
if any(~isnan(params.kiso))
  Lk(k>params.kiso)=NaN;
  Xk(k>params.kiso)=NaN;
end

% Note: should we restrict this to the upper halfplane? or will mean do
% Get the likelihood at the individual wavenumbers; average
L=-nanmean(Lk);
if isnan(L)
  % Attempt to reset
  L=1e100;
end

if nargout==3
  % Extract the moments we'll be needing for evaluation later
  df=2;
  % First should be close to df/2, second close to df/2, third is like
  % the second except formally from a distribution that should be normal
  % with mean df/2 and variance blabla/K; the last parameter is "magic"
  momx=[nanmean(Xk) nanvar(Xk) nanmean([Xk-df/2].^2)];
end

% I say, time to extract heskosl and hes here also?
  
% Get the scores at the individual wavenumbers; average
switch params.blurs
  case {0,1}
   gam=-nanmean(gammakos(k,th,params,Hk));
   % The correct gradient is too heterogeneous to be good so scale
   gam=gam.*scl;
 otherwise
  gam=NaN;
end

% Print the trajectory, seems like one element at a time gets changed
% disp(sprintf('Current theta: %8.3g %8.3g %8.3g %8.3g %8.3g',th))

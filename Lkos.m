function [Lk,Xk]=Lkos(k,th,params,Hk)
% [Lk,Xk]=Lkos(k,th,params,Hk)
%
% Computes the likelihood function for the isotropic Forsyth model in
% Olhede & Simons (2013) under UNCORRELATED loading
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The parameter vector with elements [unscaled]:
%          D   Isotropic flexural rigidity [Nm]
%          f2  The sub-surface to surface initial loading ratio
%          s2  The first Matern parameter, aka sigma^2
%          nu  The second Matern parameter
%          rho The third Matern parameter
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
% Hk       A [prod(params.NyNx)*2]-column vector of complex Fourier-domain observations
% 
% OUTPUT:
%
% Lk       A one-column vector with the wavenumbers unwrapped
% Xk       A quadratic piece of it that gets used in the analysis of residuals
%
% SEE ALSO: 
%
% LOGLIOS
%
% Last modified by fjsimons-at-alum.mit.edu, 10/22/2014

% Extract the needed parameters of the simulation variables
blurs=params.blurs;

switch blurs
 case {0,1}
  % That's lots of screen time, FMINUNC evaluates this a lot
  % disp(sprintf('%s without blurring',upper(mfilename)))
  % First calculate the Matern spectrum with the spectral parameters
  S11=maternos(k,th);

  % Then calculate then T matrices with the lithospheric parameters, and yes
  % we know Tinv will have an Inf and detT a 0 at k=0, but HFORMOS will
  % turn the Inf into a NaN and log(0)+NaN remains NaN, ...
  [invT,detT]=Tos(k,th);
  
  % Then put it all together... and all we have to worry about is a NaN in
  % Lk which we take care of in LOGLIOS. Note that Lk should be real. 
  warning off MATLAB:log:logOfZero
  Xk=hformos(S11,invT,Hk);;
  Lk=-2*log(S11)-log(detT)-Xk;
  warning on MATLAB:log:logOfZero
 otherwise
  % That's lots of screen time, FMINUNC evaluates this a lot
  % disp(sprintf('%s with blurring factor %i',upper(mfilename),blurs))

  % Blurs IS the refinement parameter; make new wavenumber grid
  [k2,kzero]=knums(params,1);
  
  % Now make the spectral-spectral portion of the spectral matrix
  S11=maternos(k2,th);
  % The lithospheric-spectral matrix on this second grid
  [~,~,~,T]=Tos(k2,th,params); 
  % Which we multiply by the spectral-spectral portion
  S=[S11.*T(:,1) S11.*T(:,2) S11.*T(:,3)];
    
  % Which we need to convolve now in two dimensions
  % And then do subsampling onto the original target grid
  Sb=bluros(S,params,1);
  
  % Now we need the determinant of the blurred S and its inverse
  detS=[Sb(:,1).*Sb(:,3)-Sb(:,2).^2];
  % If blurs then [detS-detT.*S11.^2] should be tiny
  % plot(detS,detT.*S11.^2,'+'); axis image; grid on
  % If blurs then [invS(:,1)-invT(:,1)./S11] etc should be tiny
  % plot(invS(:,1),invT(:,1)./S11,'+'); axis image; grid on
  invS=[Sb(:,3) -Sb(:,2) Sb(:,1)]./repmat(detS,1,3);
  % Trouble is at the central wave numbers, we should take those out
  
  % Then put it all together...
  warning off MATLAB:log:logOfZero
  Xk=hformos(1,invS,Hk);
  Lk=-log(detS)-Xk;
  warning on MATLAB:log:logOfZero
  
  % Behavior is rather different if this is NOT done... knowing that it
  % will not blow up but rather be some numerically large value
  Lk(kzero)=NaN;
  % Check that this is correctly done
  difer(k(kzero),[],[],NaN)
end

% Should make sure that this is real! Don't take any chances
Lk=realize(Lk);



function Sk=maternos(k,th0s,varargin)
% Sk=MATERNOS(k,[s2 nu rho],d)
%
% Calculates the three-parameter isotropic TWO-dimensional Matern
% spectral density used by Olhede & Simons (2013).
%
% INPUT:
%
% k        Wavenumber(s), e.g. from KNUM2 [rad/m] - unwrapped 
% th0s     The spectral parameter vector, containing
%          s2    The first Matern parameter, aka sigma^2
%          nu    The second Matern parameter
%          rh    The third Matern parameter
% d        Dimension [default is 2]
%
% OUTPUT:
%
% Sk       A column vector with all the wavenumbers unwrapped
%
% SEE ALSO:
%
% MATERNPRC, MATERNOS2D
%
% Last modified by fjsimons-at-alum.mit.edu, 06/27/2014

% These are always the last three elements of the input 
s2=th0s(end-2);
nu=th0s(end-1);
rho=th0s(end);

% Change dimension if you like
if nargin==3
  d=varargin{1};
else
  d=2;
end

% Adjust for dimensionality by specialization
switch d
 case 2
  pd=nu;
 otherwise
  pd=gamma(nu+d/2)/gamma(nu);
end

% Calculate the denominator in the spectral density
avark=4*nu/pi^2/rho^2+k(:).^2;
% Calculate the d-dimensional spectral density
Sk=s2*pd*nu^nu*4^nu/pi^(d/2)/(pi*rho)^(2*nu).*avark.^(-nu-d/2);

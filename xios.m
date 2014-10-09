function xi=xios(k,D,DEL,g)
% xi=xios(k,D,DEL,g)
%
% Calculates the isotropic xi in Olhede & Simons.
%
% INPUT:
%
% k          Wavenumber(s) at which this is to be evaluated [rad/m]
% D          Isotropic flexural rigidity [Nm]
% DEL        Two density contrasts, surface and subsurface [kg/m^3] 
% g          Gravitational acceleration [m/s^2]
%
% OUTPUT:
%
% phi        A column vector with all the wavenumbers unwrapped
%
% Last modified by fjsimons-at-alum.mit.edu, 04/22/2010

xi=1+D*k(:).^4/g/DEL(2);



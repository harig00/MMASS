function hfk=hformos(S,A,Hk,xver)
% hfk=hformos(S,A,Hk,xver)
%
% Computes a certain quadratic form used by Olhede & Simons, namely
% the function (1/S) [Hk^H * A * Hk] where ^H is the Hermitian transpose
%
% INPUT:
%
% S     A spectral density, which will be inverted
% A     A wavenumber-dependent SYMMETRIC matrix, e.g. Tinv, [Nx3], or
%       a wavenumber-dependent vector, [Nx1]
% Hk    A complex 2-column matrix of Fourier-domain observations, [Nx2], or
%       a complex 1-column vector of Fourier-domain observations, [Nx1]
%
% OUTPUT:
%
% hfk   A column vector with the wavenumbers unwrapped, [N]
%
% NOTE: 
%
% In our formalism, this should result in a real-valued function. 
%
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2014

defval('xver',0)

% The Inf's at k=0 get turned into NaN's here
hfk=    1./S.*[  conj(Hk(:,1)).*A(:,1).*Hk(:,1)];
if size(A,2)==3 && size(Hk,2)==2
  hfk=hfk+1./S.*[2*real(conj(Hk(:,2)).*A(:,2).*Hk(:,1))+...
		 conj(Hk(:,2)).*A(:,3).*Hk(:,2)];
end
% The 2*real is a shortcut to write the sum of the cross terms which
% should cancel any imaginary parts... it's a magnitude after all

% Still may have a tiny imaginary part... 
% Get rid of it provided it is small
hfk=realize(hfk);

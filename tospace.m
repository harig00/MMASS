function Hx=tospace(Hk,NyNx)
% Hx=TOSPACE(Hk,NyNx)
%
% Converts a vector of two-dimensional spectral-domain observations to
% the space domain with the proper normalization and again returned as a
% vector. This pertains to unitary transforms, not Matlab's standard.
%
% INPUT:
%
% Hk       A complex vector of Fourier-domain entries
% NyNx     The number of samples in the y and x directions, which equals
%          the size of the corresponding wavenumber matrix 
%
% OUTPUT:
%
% Hx       A real vector of spatial-domain entries
%
% SEE ALSO:
%
% TOSPEC
% 
% Last modified by fjsimons-at-alum.mit.edu, 01/05/2012

% This is all it is, really. Unitary transform if we
% multiply IFFT by sqrt(prod(NyNx)) and divide FFT by it. 
Hx=realize(indeks(ifft2(ifftshift(reshape(Hk,NyNx))),':')...
	   *sqrt(prod(NyNx)));

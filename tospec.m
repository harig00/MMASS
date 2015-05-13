function Hk=tospec(Hx,p)
% Hk=TOSPEC(Hx,p)
%
% Converts a vector of two-dimensional spatial-domain observations to the
% spectral domain with the proper normalization and again returned as a
% vector. This pertains to unitary transforms, not Matlab's standard. 
%
% INPUT:
%
% Hx       A real vector of spatial-domain entries
% p        A parameter structure containing
%          NyNx  The number of samples in the y and x directions which
%                equals the size of the corresponding wavenumber matrix 
%
% OUTPUT:
%
% Hk       A complex vector of Fourier-domain entries
%
% SEE ALSO:
%
% TOSPACE
% 
% Last modified by fjsimons-at-alum.mit.edu, 01/05/2012

% This is all it is, really. Unitary transform if we
% multiply IFFT by sqrt(prod(NyNx)) and divide FFT by it. 
% New thing is that we stick the actual physical dimension in here 
Hk=indeks(fftshift(fft2(reshape(Hx,p.NyNx))),':')...
   /sqrt(prod(p.NyNx))*sqrt(prod(p.dydx));

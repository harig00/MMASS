function [H,V,K,t]=choice(R,lwin)
% [H,V,K,t]=CHOICE(R,lwin)
%
% For a radius of concentration 'R' in the phase plane,
% comes up with the eigenvalues V, eigenfunctions H and
% their number K and so on. Give window length in samples.
% Eigenvalues change with R, not the functions themselves.
% The radius 'R' is in seconds knowing the time axis
% goes from -5 to +5. So in a way, the time concentration is
% 2R/10. The number of eigenvalues is ceil(R^2).
%
% Check normalization by
% diag(choice(3,100)'*choice(3,100))*10/100
% In other words, don't forget to multiply by the sampling interval, 
%
% For use in BAYRAM and COHEN

if length(lwin)~=2.^nextpow2(length(lwin)); 
error('Dyadic window length required');  end

% Take a fairly complete set of eigenvalues
K=ceil(R^2);

% Don't try to make the min and max t functions of R
% You're always multiplying by a Gaussian so go between
% its bounds.
t=linspace(-5,5,lwin);

[H,V]=hermite(0:K-1,t,R,'norm');


function detcheck(detA,A,tol)
% DETCHECK(detA,A,tol)
%
% Checks the determinant of a two-by-two SYMMETRIC matrix of wavenumbers
%
% INPUT:
%
% detA     The purported determinant [N]
% A        The actual matrix whose determinant this should be [Nx3] 
% tol      The negative of the tolerance exponent [be explicit!]
%
% Last modified by fjsimons-at-alum.mit.edu, 02/27/2012

difer(detA-[A(:,1).*A(:,3)-A(:,2).^2],tol,[],NaN)

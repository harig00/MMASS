function cholcheck(L,A,tol,meth)
% CHOLCHECK(L,A,tol,meth)
%
% Checks the Cholesky factorization of a two-by-two symmetric matrix of
% wavenumbers, noting that L is LOWER-triangular, so A=L*L^T.
%
% INPUT:
%
% L        The purported determinant Cholesky factorization [Nx3]
% A        The actual matrix whose factorization this should be [Nx3]
% tol      The negative of the tolerance exponent [be explicit~]
% meth     1 By multiplication and comparing against itself
%          2 By using Matlab's function CHOL at a random wavenumber
%
% If it passes, no output, if it doesn't, get an error
% 
% Last modified by fjsimons-at-alum.mit.edu, 02/27/2012

switch meth
 case 1
  difer(L(:,1).*L(:,1)-A(:,1),tol,[],NaN)
  difer(L(:,1).*L(:,2)-A(:,2),tol,[],NaN)
  % This element is the biggest with very little variance
  % so special rules apply to the tolerance
  difer([L(:,2).*L(:,2)+L(:,3).*L(:,3)-A(:,3)]./...
 	A(:,3)/length(A(:,3)),tol,[],NaN)
 case 2
  randi=ceil(rand*length(A)); 
  Arand=[A(randi,1) A(randi,2) ; A(randi,2) A(randi,3)];
  Lrand=[L(randi,1)    0       ; L(randi,2) L(randi,3)];
  try
    difer(chol(Arand)'-Lrand,tol,1,NaN);
  catch
    % This may have a small imaginary part
    difer(imag(Arand),tol,[],NaN)
    try 
      % But it should be positive definite if f2 is also positive  
      difer(chol(real(Arand))'-Lrand,tol,[],NaN);
    catch
      % If not, let it go but make a note of it
      disp('Cholesky not happy; most likely claims not positive definite')
    end
  end
 otherwise
  error('Specify valid method')
end

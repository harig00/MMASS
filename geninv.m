function x=geninv(A,b)
% x=GENINV(A,b)
% Ainv=GENINV(A)
%
% Gives the generalized inverse solution to Ax=b for well-behaved problems.
% b has to be a column vector.
%
% Last modified by fjsimons-at-alum.mit.edu

defval('b',[])
if isempty(b)
  x=inv(A'*A)*A';
else
  x=inv(A'*A)*A'*b;
end

% IF BADLY CONDITIONED, USE PINV


function invcheck(invA,detA,A,tol,meth,sev,goods)
% INVCHECK(invA,detA,A,tol,meth,sev,goods)
%
% Checks the inverse of a two-by-two symmetric matrix of wavenumbers 
%
% INPUT:
%
% invA     The purported inverse of A [Nx3]
% detA     The determinant of A [N]
% A        The actual matrix whose inverse this should be [Nx3]
% tol      The negative of the tolerance exponent [be explicit!]
% meth     1 By the Cayley-Hamilton theorem (uses the determinant)
%          2 By multiplication and comparing against the identity
% sev      The severity of the output message [default: 0]
% goods    Message string [default: NaN, choose [] for DIFER's default]
%
% Last modified by fjsimons-at-alum.mit.edu, 12/19/2012

% The zero wavenumbers are always going to be trouble
kbadi=find(sum(isinf(invA),2)); 
kgoodi=skip(1:size(A,1),kbadi);
% disp(sprintf('Troublesome wave number at %3g',k(kbadi)))

defval('sev',0)
defval('goods',NaN)

switch meth 
  case 1
   warning off MATLAB:divideByZero
   chek1=invA-[A(:,3) -A(:,2) A(:,1)]./repmat(detA,1,3);
   difer(chek1(~isnan(chek1))/length(chek1),tol,sev,goods)
   warning on MATLAB:divideByZero
 case 2
  chek2=invA(kgoodi,1).*A(kgoodi,1)+invA(kgoodi,2).*A(kgoodi,2)-1;
  difer(chek2(~isnan(chek2)),tol,sev,goods)
  chek3=invA(kgoodi,1).*A(kgoodi,2)+invA(kgoodi,2).*A(kgoodi,3);
  difer(chek3(~isnan(chek3)),tol,sev,goods)
  chek4=invA(kgoodi,2).*A(kgoodi,2)+invA(kgoodi,3).*A(kgoodi,3)-1;
  difer(chek4(~isnan(chek4)),tol,sev,goods)
 otherwise
  error('Specify valid method')
end


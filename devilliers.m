function [D,V,g,c,N,A,D2,V2,g2]=devilliers(m,N,K,L)
% [D,V,g,c,N,A,D2,V2,g2]=DEVILLIERS(m,N,K,L)
%
% Builds the tridiagonal matrix that we need to diagonalize in order to
% find the Jacobi/Bessel expansion coefficients of the radial part of the
% Slepian functions concentrating to a disk-shaped region.
%
% INPUT:
%
% m    The fixed order of the concentration eigenfunctions [default: 3]
% N    The Shannon number of the full concentration problem [default: 3]
% K    The number of functions required [default: N]
% L    The maximum degree of the Jacobi/Bessel expansion [default: 2N]
%
% OUTPUT:
%
% D    The expansion coefficients, which is an (L+1)x(L+1) matrix with
%      the row index the degree and column the eigenfunction rank 
% g    The integral equation eigenvalues
% V    The concentration eigenvalues
% c    The concentration factor c=2*sqrt(N)
% N    The Shannon number of the concentration problem
% A    The tridiagonal matrix whose eigenvectors are D
% D2   An alternative orthogonal set of coefficients from a symmetric A
% V2   An alternative set of concentration eigenvalues belonging to D2
% g2   An alternative set of integral equation eigenvalues
%
% SEE ALSO: SWDISK
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/04/2013
% Last modified by dongwang-at-princeton.edu, 08/06/2008

% This is probably equivalent to the three-term recursion of Slepian
% (1964). See also Shkolnisky (2007) and Xiao (2001) for the m=1/2
% case, i.e. for the "regular" one-dimensional prolate spheroidals.

% Supply defaults just in case, though usually this is a called routine 
defval('m',3)
defval('N',3)
defval('K',ceil(N))
% Some heuristic defaults - make identical in SWDISK
defval('L',min(max([2*K ceil(2*N) 2*m 10]),84))

% Stick to the notation of Slepian (1964)
c=2*sqrt(N);
disp(sprintf('c = %4.2g ; m = %i ; N = %4.2g ; K = %i ;  L = %i ',...
	     c,m,N,K,L))

% Where does the range start? From zero.
k=0:L; kk=0:L-1;

% Make the sub, main and superdiagonal
warning off MATLAB:divideByZero
ldiag=-c^2*(m+kk+1).^2./(2*kk+m+1)./(2*kk+m+2);
ondiag=(2*k+m+1/2).*(2*k+m+3/2)+...
       c^2/2*(1+m^2./(2*k+m)./(2*k+m+2));
udiag=-c^2*(  kk+1).^2./(2*kk+m+2)./(2*kk+m+3);
warning on MATLAB:divideByZero

% Form the tridiagonal non-symmetric matrix
A=tridiag(ldiag,ondiag,udiag);

% Now for the symmetric variety. Common form for both diagonals
uldiag=-c^2*(kk+1).*(m+kk+1)./sqrt(2*kk+m+1)./(2*kk+m+2)./sqrt(2*kk+m+3);
% And notice that A and AA and AAA have the same eigenvalues! But now, how do the
% eigenvectors scale? 
AA=tridiag(uldiag,ondiag,uldiag);

% Supply the correct entry when the order was zero
if m==0
  A(1,1)=1/2*3/2+c^2/2;
  AA(1,1)=A(1,1);
end

% Now find the eigenvectors of this thing
if K>=L & K<=L+1
  % Find them all - somewhat slower
  [D,v]=eig(A);
  % But when you solve the symmetric system you get orthogonal eigenvectors!
  [D2,v2]=eig(AA);
  difer(sort(v(:))-sort(v2(:)),[],[],NaN)
elseif K<L
  % Return the smallest eigenvalues of the Sturm-Liouville system; these
  % correspond to the largest eigenvalues of the concentration problem
  OPTS.disp=0;
  [D,v]=eigs(A,K,'SM',OPTS);
  [D2,v2]=eigs(AA,K,'SM',OPTS);
  difer(sort(v(:))-sort(v2(:)),[],[],NaN)
elseif K>L+1
  error(sprintf('You must increase truncation degree to at least %i',K))
end

% Note that the columns of D are the eigenvectors; that they are not
% necessarily an orthogonal set; that each column has been normalized to
% energy one; so when using them in the expansion, we will still need to
% normalize the functions according to our particular wishes of the day.

% Note that the sorting may be non-existent or wrong
[v,i]=sort(diag(v)); D=D(:,i);

% The rows of D are the degrees l, the columns the rank of the function
% Now find the concentration eigenvalues as per Slepian (1964) eq. 46
g=c^(m+1/2)*D(1,1:K)./2^(m+1)./factorial(m+1)./sum(D(:,1:K),1);
V=g.^2*c;
g2=c^(m+1/2)*D2(1,1:K)./2^(m+1)./factorial(m+1)./sum(D2(:,1:K),1);
V2=g2.^2*c;

% Check that the eigenvalues are properly contained
if ~isnan(nanmean([V(V>1) V(V<0)]));
  V1=V-1; V1=V1(V1>0); V0=V+1; V0=V0(V0<1);
  if nanmean([V1 V0])>1000*eps
    error('Eigenvalues exceeding 0 to 1 - definitely increase L'); 
  end
end

% Build in a check on the convergence - we want the last term in the row
% sum to not matter at all for the first K requested functions (noting
% that technically, we may have requested K=L+1).
avco=abs(D(L+1,1:K))./sum(abs(D(1:L,1:K)),1);

% The tough criterion is when the sum converges, period
tough=mean(abs(avco(1:K)));
%disp(sprintf('Tougher convergence estimate %6.3e',tough))

% The lenient criterion is when the sum converges for the big eigenvalues
lenient=mean(abs(avco(g(1:K)>eps)));
%disp(sprintf('Lenient convergence estimate %6.3e',lenient))

% Formal error reporting on these two criteria
difer(tough,[],[])
difer(lenient,[],[])


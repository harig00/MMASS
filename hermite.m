function [H,V]=hermite(nvec,x,varargin)
% H=hermite(n,x)
% H=hermite(n,x,[],'norm')
% [H,V]=hermite(n,x,R)
% [H,V]=hermite(n,x,R,'norm')
%
% Computes (normalized) Hermite polynomials with compact support.
% H Hermite eigenvectors and V their eigenvalues of a disk-shaped domain
% in the time-frequency phase plane of radius R.
%
% Check normalization by
% diag(choice(3,100)'*choice(3,100))*10/100
% In other words, don't forget to multiply by the sampling interval, 
%
% If n is a vector returns polynomials in the columns of H.
%
% See also LAGUERRE, COHEN
%
% Example I 
% 
% x=linspace(-5,5,100); plot(x,hermite(0:4,x,[],'norm'),'k-')
% H=hermite(0:4,x,[],'norm'); 
% H'*H*indeks(diff(x),1) % shows these are orthonormal
% 
% Example II
% 
% x=linspace(-2,2,100); plot(x,hermite(0:4,x),'k-')
%
% Example III
%
% n=0:20; [H,V]=hermite(n,x,3,'norm'); plot(n,V,'o-'); hold on
% [H,V]=hermite(n,x,5,'norm'); plot(n,V,'o-');
% [H,V]=hermite(n,x,7,'norm'); plot(n,V,'o-');

% $H_n(x)$; or the orthonormal
% $\pi^{-1/4}(2^nn!)^{-1/2}e^{-x^2/2}H_n(k)$ as in Bayram and Baraniuk (2000)
%
% http://mathworld.wolfram.com/HermitePolynomial.html
%
% \begin{eqnarray} 
% H_n(x) &=& (-1)^ne^{x^2} {d^n\over dx^n} e^{-x^2}\\
%            &=& e^{x^2/2}\left({x - {d\over dx}}\right)^n e^{-x^2/2}
% \end{eqnarray}
%
%
% Using recurrence relation
% $H_{n+1}(x)=2xH_n(x)-2nH_{n-1}(x)$

% Daubechies, IEEE Inf. Theor., 34 (4) 605-612, 1988
% Bayram and Baraniuk, Proc. IEEE Sig. Proc. 173-176, 1996

% Last modified by  FJS, Aug 31th 2000

nvec=nvec(:);

for index=1:size(nvec)
  n=nvec(index);
  switch n
    case 0
      H=ones(size(x));
    case 1
      H=2*x;
    otherwise
      Hnm2=1;  % H_0 
      Hnm1=2*x; % H_1
      for m=1:n-1            
	H=2*x.*Hnm1-2*m*Hnm2; % This is now $H_{m+1}$
	Hnm2=Hnm1;
	Hnm1=H;
      end
  end
  % Normalization
  if nargin==4 & varargin{2} == 'norm'
    gauz=exp(-x.^2/2);
    H=H*pi^(-1/4)/sqrt(2^n*gamma(n+1)).*gauz;
  end
  Hout(:,index)=H(:);
end

H=Hout;

if nargout==2 & ~isempty(varargin{1})
  V=gammainc(varargin{1}^2/2,nvec+1);
end



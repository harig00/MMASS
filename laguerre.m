function L=laguerre(nvec,k,x,varargin)
% L=LAGUERRE(n,k,x)
% L=LAGUERRE(n,k,t,f)
% L=LAGUERRE(n,k,x,[],'bayram')
% L=LAGUERRE(n,k,t,f,'bayram')
%
% Computes associated Laguerre polynomials 
% 'n' can be a vector; 'k=0' gives unassociated polynomials 
%
% See HERMITE
%
% Example I:
%
% t=linspace(-1,5,100);
% L=laguerre(0:5,0,t,[],'bayram'); plot(t,L,'-k')
% L'*L*indeks(diff(t),1) % These are not orthonormal
%
% f=linspace(-1,1,100);
% imagesc(sum(L,3)); axis image
%
% Example II:
%
% plot(linspace(-1,5,100),laguerre(0:5,0,linspace(-1,5,100)))

% $L_n^k(x)$ for a 1-D domain;
% $L_n^k(t^2+f^2)$ on a 2-D domain; or
% $e^{-\pi/2[t^2+f^2]}L_n^k(\pi[t^2+f^2])$ as in Bayram and Baraniuk (2000)
%
%
% http://mathworld.wolfram.com/LaguerrePolynomial.html
%
% \begin{eqnarray} 
% L_n^k(x) &=& {e^xx^{-k}\over n!}{d^n\over dx^n}(e^{-x}x^{n+k})\\
%          &=& \sum_{m=0}^n (-1)^m {(n+k)!\over (n-m)!(k+m)!m!}x^m, 
% \end{eqnarray}

% Last modified by  FJS, Aug 31th 2000

if nargin>=4 & ~isempty(varargin{1})
  M=length(x); N=length(varargin{1});
  [T,F]=meshgrid(x,varargin{1});
  x=T.^2+F.^2;
end  

x=x(:)';
for index=1:length(nvec)
  n=nvec(index);
  m=0:n;
  lb=gamma(n+k+1)./(gamma(k+m+1))./gamma(n-m+1)./gamma(m+1).*(-1).^m;
  if nargin==5 & varargin{2}== 'bayram'
    rb=[ones(size(x)) ; cumprod(repmat(pi*x,n,1),1)];
    L(:,index)=[exp(-pi/2*x).*(lb*rb)]';
  else
    rb=[ones(size(x)) ; cumprod(repmat(x,n,1),1)];
    L(:,index)=[lb*rb]';  
  end
end

if nargin>=4 & ~isempty(varargin{1})
  for index=1:length(nvec)
    Lout(:,:,index)=reshape(L(:,index),M,N);
  end
  L=Lout;
end


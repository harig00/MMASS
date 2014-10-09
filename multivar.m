function [Sabll,K,V]=multivar(L,TH,sord,l)
% [Sabll,K,V]=MULTIVAR(L,TH,sord,l)
%
% Calculates the single-degree multitaper covariance matrix.
% Dahlen & Simons (2008) eq. (163)
%
% INPUT:
%
% L         Taper bandwidth that is required in this case [scalar]
% TH        Colatitudinal radius of the cap, in degrees <=180 [scalar]
%           Colatitudinal halfwidth of the cut, degrees <=90 [scalar]
% sord      1 Single cap of diameter 2TH [default]
%           2 Double cap left by subtracting belt of width 2TH
% l         The degree that you are looking at [scalar; default: 2*L]
%
% Calculates multitaper covariance matrix: at a single degree, but for
% all the relevant multitapers arranged in a matrix of decreasing
% concentration eigenvalue. This is in the white signal and white noise,
% or moderately colored approximation, all signal and noise strength
% factors are set to or sum to one.
%
% OUTPUT:
%
% Sabll     The covariance between differently tapered estimates at l
% K         The Shannon number of the problem
% V         The eigenvalues of the concentration problem
%
% SEE ALSO: PERIODOVAR, GAMMAB, GAMMAP, WHITEVARIANCE, MTVAR
%
% Last modified by fjsimons-at-alum.mit.edu, 04/25/2007

defval('L',15)
defval('TH',30)
defval('sord',1)
defval('l',2*L)

if length(l)>1
  error('Degree l must be a scalar')
end

Lpot=(L+1)^2;

% First get the handy GAMMAB matrix at the even degrees only
[Gabp,p,K,V]=gammab(L,TH,sord,1);

disp(sprintf('We are getting a rounded Shannon number of %i',round(K)))

% Make this properly (L+1)^2 x (L+1)^2 x (L+1)
Gabp=squeeze(Gabp);

% Better sort them in function of decreasing eigenvalue
[V,j]=sort(V,'descend');

% This is now properly sorted
Gabp=Gabp(j,j,:);

% Get the Wigner 3j symbols with bottom-row of zero
% Perhaps later change this to use ZEROJ instead
W=wigner0j(2*L,l,l);

% Only select the evens since we're doing variance at equal l=l'
W2=(2*p+1).*W(1:2:end).^2;

% Note every layer of Gabp is symmetric in ab
Sabll=zeros(Lpot,Lpot);
for indr=1:Lpot
  for indc=indr:Lpot
    Sabll(indr,indc)=W2*squeeze(Gabp(indr,indc,:))/2/pi;
  end
end

% Symmetrize at the end in one step
Sabll=Sabll+tril(Sabll',-1);

% Best is to verify this with what MVARRATIOS returned! [MT variance]
Smtll=V*Sabll*V'/K^2;

% Whole-sphere variance! [WS variance]
Swsll=2/(2*l+1);

% Check this is what we have been thinking it was
if L<=20
  difer(Smtll/Swsll-mvarratios(L,TH,sord,l))
end


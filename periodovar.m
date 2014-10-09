function [v1,v2,v3,v4,TH,A,Bl,bels]=periodovar(l,TH,Bl,bels,sord,meth,Lmax)
% [v1,v2,v3,v4,TH,A,Bl,bels]=PERIODOVAR(l,TH,Bl,bels,sord,meth,Lmax)
%
% Calculates the ratio of the variance of the periodogram estimate to the
% whole-sky estimate; exact for white spectra and noise.
% Dahlen and Simons (2008) eq. (174-175)
%
% INPUT:
%
% l         The spherical harmonic degree where this is valid [scalar]
% TH        Size of the polar cap, in degrees; may be vector
% Bl        The power spectrum of the polar caps, if available
% bels      The degrees at which this power spectrum is evaluated
% sord      1 Single cap of diameter 2TH [default]
%           2 Double cap left by subtracting belt of width 2TH
%           3 Equatorial belt of width 2TH
%
% --------> Options pertaining to the calculation of needed power spectra:
% meth      1 Using unnormalized Legendre polynomials [default]
%           2 Integrating normalized Legendre polynomials
% Lmax      Maximum degree to which you'll be taking this... [default: 2*l]
%
% OUTPUT:
%
% v1         Ratio of the variances as a function of degree
% v2         An asymptotic approximation for large l, independent of l
% v3         An asymptotic approximation for small TH, independent of l
% v4         An asymptotic approximation for large TH, independent of l
% TH         The actual polar cap sizes used
% A          The area of the actual polar caps
% Bl         The power spectrum of the polar caps
% bels       The degrees at which this power spectrum is evaluated
%
% SEE ALSO: MULTIVAR
%
% Last modified by fjsimons-at-alum.mit.edu, 04/29/2007

defval('l',60)
defval('TH',[0.1 2:2:180])
defval('xver',0)
defval('Bl',[])
defval('sord',1)
defval('meth',1)

if length(l)>1
  error('Degree l must be a scalar')
end

% It's the variance, only the even degrees of the boxcar ever play a role
evens=1;
% And we knew the Lmax we need to succeed...
% Only if we want to return Bl and reuse it do we need to specify this
defval('Lmax',2*l);

% For the areas, see also BCOUPLING - or maybe make new function altogether?
A=4*pi*spharea(TH,sord);

if isempty(Bl)
  % Get the power spectrum of the boxcar taper
  [Bl,bels]=bpboxcap(TH,Lmax,meth,evens,sord);
end

% Get the Wigner 3j symbols with bottom-row of zero
W=wigner0j(Lmax,l,l);

% If you want to be excessive about verifying this
if xver==1
  difer(W-wigner3jm(Lmax,l,l,0,0,0));
  % If more bandwidth database should do threej too
  disp('PERIODOVAR excessive verification passed')
end

% Only select the evens since we're doing variance at equal l=l'
W=W(1:evens+1:end);

% Put them all together and factor in the area
v1=sum(Bl.*repmat(sqrt(2*bels'+1).*W(:),1,length(TH)).^2);
v1=v1*(2*l+1)*4*pi./A.^2;
% Check this for both methods of power spectral computations
% For the tiniest areas only this may make a difference

% Now make an approximation valid for large l - the curves tend to this
% value - the xlm respect whatever choice you'd made about evenness
v2=(4*pi)^2./A.^2.*sum(Bl.*(xlm(bels,0,repmat(pi/2,1,length(TH))).^2));

% THESE here aren't so good
% ... and small TH (note that then A=pi*TH^2 for single caps only)
v3=32/3/pi^2*sqrt(pi./A);
% ... and large TH
v4=1+2/3/pi*(pi-TH*pi/180).^3;

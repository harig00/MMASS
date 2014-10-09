function varargout=transl(f,Te,r,dr,E,v,g)
% [k12,l12,E,v]=TRANSL(f,Te,r,dr,E,v,g)
%
% Calculates the coherence-square halfway point in function 
% of Forsyth's subsurface-to-surface loading ratio and the elastic
% thickness in km. Returns the wavenumber at which C2=0.5, in rad/km.
%
% INPUT:
%
% f         Amplitude loading ratio [dimensionless, not squared]
% Te        Elastic thickness [km], can be vector
% r         Density [kg/m3]
% dr        Density contrast [kg/m3]
% E         Young's modulus [Pa]
% v         Poisson's ratio [dimensionless]
% g         Gravitational acceleration [m/s2]
%
% OUTPUT:
%
% k12       Wavenumber at which C2=0.5, in rad/km
% l12       Wavelength at which C2=0.5, in km
% E         Young's modulus
% v         Poisson's ratio
%
% SEE ALSO: FORSYTH, COHERETE, COHEREF2
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/05/2012

defval('f',0.5)
defval('Te',17)
defval('E',1.4e11);
defval('v',0.25);
defval('r',2670);
defval('dr',630);
defval('g',fralmanac('GravAcc'));

% Convert to standard units [m]
Te=Te*1000;

disp(sprintf('E= %5.3g; v= %5.3f',E,v))

% How to get to this point
% coh=[ '((1+D*k^4/dr/g)*dr^2+f^2*r^2*(1+D*k^4/r/g))^2/...
%        ((1+D*k^4/dr/g)^2*dr^2+f^2*r^2)/...
%        (dr^2+f^2*r^2*(1+D*k^4/r/g)^2)'];
% ss=[ 'solve(0.5=' coh ',k);'];
% maplesol=parse(maple(ss),',');
% k12=eval(maplesol(1,:))

D=(E*Te.^3)./(12*(1-v^2));

[FF,DD]=meshgrid(f,D);
k12=1./2.*2.^(3./4)./DD./FF.*...
    ((-FF.*dr-FF.*r+FF.^2.*r+dr+...
      (FF.^2.*dr.^2+4.*FF.^2.*dr.*r-2.*FF.^3.*dr.*r+...
       2.*FF.*dr.^2+FF.^2.*r.^2+2.*FF.^3.*r.^2-...
       2.*FF.*r.*dr+FF.^4.*r.^2+dr.^2).^(1./2))...
     .*g.*DD.^3.*FF.^3).^(1./4);

% Now do it exactly as it appears in the EPSL paper:
beta=(FF.^2.*dr.^2+4.*FF.^2.*dr.*r-2.*FF.^3.*dr.*r+...
       2.*FF.*dr.^2+FF.^2.*r.^2+2.*FF.^3.*r.^2-...
       2.*FF.*r.*dr+FF.^4.*r.^2+dr.^2);

% Next line has the typo
%k12=(2*g./DD./FF.*(-FF.*dr-FF.*r+FF.^2.*r+dr+beta.^(1/2))).^(1./4);

% And noting the typo, this should be corrected to
k=(1/2*g./DD./FF.*(-FF.*dr-FF.*r+FF.^2.*r+dr+beta.^(1/2))).^(1./4);

% Check the results
difer(k-k12,[],[],NaN)

% Also check by plugging it in
difer(((1+DD.*k.^4/dr/g)*dr^2+FF.^2*r^2.*(1+DD.*k.^4/r/g)).^2./...
       ((1+DD.*k.^4/dr/g).^2*dr^2+FF.^2*r^2)./...
        (dr^2+FF.^2*r^2.*(1+DD.*k.^4/r/g).^2)-1/2,[],[],NaN)

% Back into the right units
k12=k12*1000;
l12=2*pi./k12;

% Produce output
varns={k12,l12,E,v};
varargout=varns(1:nargout);

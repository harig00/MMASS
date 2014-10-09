function I=sinsinexp(m1,m2,mp)
% I=SINSINEXP(m1,m2,mp)
%
% INPUT: 
%
% m1,m2,m3     Three spherical harmonic orders
% 
% Calculates the integral of 
%
% \int_{0}^{pi}\sin(m1\phi)\sin(m2\phi) exp(i*mp*phi)\,d\phi
%
% Silly me - this integral is no good: usually, we want 2*pi
%
% Last modified by fjsimons-at-alum.mit.edu, 08.08.2006

cospim1=cos(pi*m1); cospim2=cos(pi*m2);
sinpim1=sin(pi*m1); sinpim2=sin(pi*m2);
expimp=exp(mp*pi*i);

itop=-[-expimp*mp*m1^2*sinpim1*sinpim2*i-expimp*mp*m2^2*sinpim1*sinpim2* ...
       i-2*i*expimp*mp*m1*m2*cospim1*cospim2+expimp*mp^3*sinpim1*sinpim2* ...
       i-expimp*m1*m2^2*cospim1*sinpim2-expimp*m2*mp^2*sinpim1*cospim2- ...
       expimp*mp^2*m1*cospim1*sinpim2-expimp*m1^2*m2*sinpim1*cospim2+ ...
       expimp*m2^3*sinpim1*cospim2+expimp*m1^3*cospim1*sinpim2+2*i*mp*m1* ...
       m2];
ibot=[mp^4-2*mp^2*m1^2-2*mp^2*m2^2+m1^4-2*m1^2*m2^2+m2^4]; 

% Limiting case Using L'Hopital
if abs(ibot)<eps
  disp('Using L''Hopital on m1')
  % Derivative wrt m1; could use m2 or mp, why not?
  itop=[-2*i*mp*m2+expimp*mp*m1^2*cospim1*pi*sinpim2*i-expimp*mp^3*cospim1* ...
	pi*sinpim2*i-2*i*expimp*mp*m1*m2*sinpim1*pi*cospim2+2*i*expimp* ...
	mp*m1*sinpim1*sinpim2+expimp*mp*m2^2*cospim1*pi*sinpim2*i+expimp* ...
	m2^2*cospim1*sinpim2-expimp*m1*m2^2*sinpim1*pi*sinpim2+expimp*m2*mp^2* ...
	cospim1*pi*cospim2+expimp*mp^2*cospim1*sinpim2-expimp*mp^2*m1* ...
	sinpim1*pi*sinpim2+2*expimp*m1*m2*sinpim1*cospim2+expimp*m1^2*m2* ...
	cospim1*pi*cospim2-expimp*m2^3*cospim1*pi*cospim2-3*expimp*m1^2* ...
	cospim1*sinpim2+expimp*m1^3*sinpim1*pi*sinpim2+2*i*expimp*mp*m2* ...
	cospim1*cospim2];
  ibot=[-4*mp^2*m1+4*m1^3-4*m1*m2^2];
end

I=itop/ibot;






function [J,coordd,dxi,deta]=cubejac(neta,nxi,sc,np)
% [J,coordd,dxi,deta]=cubejac(neta,nxi,sc,np)
% 
% Calculates the Jacobian and the coordinates of the grid nodes, not the
% pixel centers, of a single face of a cubed sphere a la Ronchi et al.,
% doi:10.1006/jcph.1996.0047 and further modified by Loris et al.
%
% INPUT:
%
% nxi,neta     Number of elements in the xi,eta direction
% sc           0 regular cubed sphere [default]
%              1 superchunk cubed sphere
% np           1 Grid node registration [default]
%              2 Pixel center registration [not ready yet]
%
% OUTPUT:
%
% J            The Jacobian of this transformation
% coordd       The Cartesian coordinates of a single face
% dxi,deta     The angular spacings
%
% SEE ALSO:
%
% CUBE2SPHERE, ANGULARPART513TO128
%
% Last modified by fjsimons-at-alum.mit.edu, 01/20/2011

% First make grid vectors of the angular part
defval('neta',2^6+1)
defval('nxi',2^6+1)
defval('sc',0)
defval('np',1)

% First make grid vectors on the nodes
if sc==0
  edgy=[-pi/4 pi/4]; 
elseif sc==1
  edgy=[-3*pi/8 3*pi/8]; 
end
  
rngy=edgy(2)-edgy(1);
xi =linspace(edgy(1),edgy(2),nxi);
eta=linspace(edgy(1),edgy(2),neta);
dxi=rngy/(nxi-1);
deta=rngy/(neta-1);

% Then make the actual grid
[ETA,XI]=meshgrid(eta,xi);

% Then make the Jacobian and the coordinates
denom=sqrt(1+tan(XI).^2+tan(ETA).^2);

% Make the Jacobian - this is \sin(\theta)\partial\theta/\partial\eta
J=(1+tan(XI).^2).*(1+tan(ETA).^2)./denom.^3;

% Make the Cartesian coordinates of the face that will be rotated
coordd=[1./denom(:) tan(XI(:))./denom(:) tan(ETA(:))./denom(:)]'; 
%coordd=[1./denom(:) tan(ETA(:))./denom(:) -tan(XI(:))./denom(:)]'; 


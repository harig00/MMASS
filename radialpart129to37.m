function xp=radialpart129to37(x)
% xp=RADIALPART129TO37(x)
%
% INPUT:
%
% x         The three-dimensional array of size MxNx129
%
% OUTPUT:
%
% xp        The three-dimensional array of size MxNx37 which also
%           includes the Jacobian already
%
% Assuming 129 equally spaced radii (edges), numerically integrate to 37
% layers, according to hardwired layers provided in here.

% Written by Ignace Loris (igloris@vub.ac.be) on 19/06/2009
% Last modified by fjsimons-at-alum.mit.edu, 07/14/2009

% Error check on hardwired dimensions
nx3=size(x,3);
if nx3~=129
  error('Input should have size 129 in third dimension!')
end

% Prepare array
xp=zeros(size(x,1),size(x,2),37);

% Hardwire the radii, in km, and produce the jacobian, in km^2
% The Earth's total radius
R_E=6371;
% The Earth's core radius (the core-mantle-boundary)
R_CMB=3481.4;
% The thickness of the layers
dr=(R_E-R_CMB)/(nx3-1);
% The location of the interfaces themselves
r=linspace(R_CMB,R_E,nx3);
% Produce the Jacobian and include the sampling interval
Jdr=r.^2*dr;
% Multiply this into the input matrix already
prefactor=repmat(Jdr,[size(x,2) 1 size(x,1)]);
% Special provision if x is effecively unidimensionsl
if size(x,1)~=1 && size(x,2)~=1
  prefactor=shiftdim(prefactor,2);
else
  prefactor=reshape(prefactor,size(x));
end
x=x.*prefactor;

% Perform the numerical integration:
% Use trapezoid rule for 22.575 km layers:
xp(:,:,37)=(x(:,:,129)+x(:,:,128))/2;
xp(:,:,36)=(x(:,:,128)+x(:,:,127))/2;

% Use Simpson's rule for 45.15 km layers:
xp(:,:,35)=(x(:,:,125)+4*x(:,:,126)+x(:,:,127))/3;
xp(:,:,31)=(x(:,:,111)+4*x(:,:,112)+x(:,:,113))/3;
xp(:,:,30)=(x(:,:,109)+4*x(:,:,110)+x(:,:,111))/3;
xp(:,:,27)=(x(:,:, 99)+4*x(:,:,100)+x(:,:,101))/3;
xp(:,:,26)=(x(:,:, 97)+4*x(:,:, 98)+x(:,:, 99))/3;
xp(:,:, 2)=(x(:,:,  3)+4*x(:,:,  4)+x(:,:,  5))/3;
xp(:,:, 1)=(x(:,:,  1)+4*x(:,:,  2)+x(:,:,  3))/3;

% These layers are the thickest: 90.3 km
% Use extended Simpson's rule, with weights [1 4 2 4 1]/3;
for i=2:24
  xp(:,:,i+1)=(x(:,:,4*i-3)+4*x(:,:,4*i-2)+2*x(:,:,4*i-1)+4*x(:,:,4*i)+x(:,:,4*i+1))/3;
end
for i=26:27
  xp(:,:,i+2)=(x(:,:,4*i-3)+4*x(:,:,4*i-2)+2*x(:,:,4*i-1)+4*x(:,:,4*i)+x(:,:,4*i+1))/3;
end
for i=29:31
  xp(:,:,i+3)=(x(:,:,4*i-3)+4*x(:,:,4*i-2)+2*x(:,:,4*i-1)+4*x(:,:,4*i)+x(:,:,4*i+1))/3;
end


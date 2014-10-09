function Z=planef(abc,x,y,plotornot)
% Z=PLANEF(abc,x,y,plotornot)
%
% Draws the plane defined by 'Z'=a*X+b*Y+c on the domain defined
% by X and Y, which are the meshgridded input 'x' and 'y'
% Written for 'planefit'. Plots or not according to 'plotornot'
%
% Last modified by fjsimons-at-alum.mit.edu, September 26th 1998

% Make sure x and y are row vectors
x=x(:)';
y=y(:)';
[a,b,c]=deal(abc(1),abc(2),abc(3));

[X,Y]=meshgrid(x,y);

Z=a*X+b*Y+c;

if plotornot
  surf(X,Y,Z)
end


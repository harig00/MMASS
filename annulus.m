function lodge=annulus(matrix,r,dr,theta,dtheta,senter)
% lodge=ANNULUS(matrix,r,dr,theta,dtheta)
% lodge=ANNULUS(matrix,r,dr,theta,dtheta,senter)
%
% Finds the coordinate indices of 'matrix' of those elements contained
% in  a band given by '(r,theta)' and '(dr,dtheta)'
% Note: we're using atan2 so the results are within -pi and pi.
% So the complete circle is defined as from th to th+dth or (-pi,2*pi)
%
% Option to define something else but the geometric center [cx cy]
%
% Example:
%
% mat=peaks(300); 
% mat(annulus(mat,120,10,pi/3,pi/6))=NaN;
% imagesc(mat) ; axis image
%
% lodge=annulus(mat,80,20,-2*pi/3,2*pi/6);
% mat(lodge)=mean(mat(lodge));
% imagesc(mat) ; axis image

% Last modified by fjsimons-at-mit.edu, May 22nd, 2001

[m,n]=size(matrix);

if nargin==5 | isempty(senter)
  cy=ceil(m/2);
  cx=ceil(n/2);
else
  cx=senter(1);
  cy=senter(2);
end

[X,Y]=meshgrid(1:n,1:m);
X=X-cx;
Y=cy-Y;

R=sqrt(X.^2+Y.^2);
TH=atan2(Y,X);

lodge=(R>=r & R< r+dr & TH>theta & TH<=theta+dtheta);

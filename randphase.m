function rp=randphase(varargin)
% This function returns random phases 
% (distributed over the interval -pi to pi)
%
% function rp=RANDPHASE(m,n)
%             (returns a m by n matrix)
%          rp=RANDPHASE(n)
%             (returns a column vector of length n)
%          rp=RANDPHASE
%             (returns one number)
%

% After Lowry

% Written by FJS, December 20th 1998

% Random numbers between 0 and 1
switch nargin
  case 2
    rp=rand(varargin{1},varargin{2});
  case 1
    rp=rand(varargin{1},1);    
  case 0
    rp=rand(1);
end

isn=rp>0.5;

% Random numbers between 0 and 0.5
rp=(-1).^isn.*rp+isn*1;

% Polynomial coefficients (from high to low)
C1=[0.010328 0.802853 2.515517];
C2=[0.001308 0.189269 1.432788 1];

% Random numbers between -pi and pi
t=sqrt(log(1./rp.^2));
rp=(-1).^isn.*(t-polyval(C1,t)./polyval(C2,t));





function varargout=Lsurfos(ND,Nf,th,params,Hk,k)
% [L2D,DD,ff2]=LSURFOS(ND,Nf,th,params,Hk,k)
%
% Computes the likelihood surface around a certain target point
%
% INPUT:
%
% ND       The size of the parameter "grid" for D
% Nf       The size of the parameter "grid" fo f2
% th       The five-parameter vector argument
%          th(1)=D    Isotropic flexural rigidity [Nm]
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% Hk       A complex matrix of Fourier-domain observations
% k        the wavenumbers at which these are being evaluated [1/m]
%
% OUTPUT:
%
% L2D      The two-dimensional likelihood around D,f2
% DD       The grid in D
% ff2      The grid in f2
%
% EXAMPLE:
%
% Lsurfos('demo1')
%
% Last modified by fjsimons-at-alum.mit.edu, 12/20/2010

defval('ND',10)
if ~isstr(ND)
  % Expand the vector for passing onto the single-wavenumber ones
  D=th(1);
  f2=th(2);
  
  % Define a grid on which to work
  rng=25;
  DD=logspace(log10(D/rng),log10(rng*D),ND);
  ff2=linspace(0,10*f2,Nf);
  
  % Calculate the likelihood surface... intensively, at first
  L2D=nan(ND,Nf);
  for ind=1:ND
    for jnd=1:Nf
      L2D(ind,jnd)=-nanmean(Lkos(k,[DD(ind) ff2(jnd) th(3:5)],params,Hk));
    end
  end

  % Output
  varns={L2D,DD,ff2};
  varargout=varns(1:nargout);
elseif strcmp(ND,'demo1')
  ND=100;  Nf=100;
  [Hx,Gx,th0,k,Hk,Gk,dydx,NyNx,DEL,g,z2]=simulos;
  [L2D,DD,ff2]=Lsurfos(ND,Nf,th0,params,Hk,k);
  % Make the figure
  clf
  imagesc(DD,ff2,L2D)
  hold on
  pth=plot(th0(1),th0(2),'o');
  hold off
  set(pth,'MarkerF','w')
end

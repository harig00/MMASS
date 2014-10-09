function varargout=forsyth(Te,lbd,f2,r,rc,drho,T,g,xver)
% [G2b,k,l,Zb,Zf]=FORSYTH(Te,lbd,f2,r,rc,drho,T,g,xver)
%
% Calculates the predicted coherence-square and admittance between
% Bouguer/free air anomalies and topography over a range of wavenumbers for a
% lithosphere loaded at the surface and at one other depth which is also
% the depth of compensation. For the coherence it doesn't matter where
% this depth is; it is implicitly included in the loading factor.
% See Forsyth (1985).
%
% Works for a vector of Te OR a vector of f2 OR a vector of r (only one vector).
%
% INPUT:
%
% Te      elastic thickness of plate [m]
% lbd     the wavelengths considered [m] [defaulted to reasonable values]
% f2      spectral ratio of bottom-to-top applied loads [dimensionless, squared]
%         0 for surface loading only (produces unity coherence)
%         1 for equal loads
%         infinity for Moho loading only (produces unity coherence)
% r       initial-load correlation coefficient
% rc      density contrast across topography interface [kg/m^3]
%         rho_crust for loading of continents
% drho    density contrast at compensation interface
%         (usually = rho_mantle - rho_crust) [kg/m^3]
% T       depth to the density contrast/compensation interface [m]
% g       Gravitational acceleration (in m/s^2) [defaulted]
%
% OUTPUT:
%
% G2b       Bouguer coherence-square (between 0 and 1)
% k         wavenumber [rad/m]
% l         wavelength [m]
% Zb        Bouguer addmittance, in mgal/m
% Zf        Free-air addmittance, in mgal/m
%
% See also MCKENZIE, ADMITTANCE, TRANSL
%
% EXAMPLE:
%
% forsyth('demo')
%
% Last modified by fjsimons-at-alum.mit.edu, 03/07/2012

defval('Te','demo')

if ~isstr(Te)
  defval('Te',[20 80]*1e3);
  defval('lbd',linspace(10e3,2000e3,100))
  defval('f2',1);
  defval('rc',2670);
  defval('drho',630);
  defval('T',35e3)
  % Young's modulus
  defval('E',1.4e11);
  % Poisson's ratio
  defval('v',0.25);
  % Gravity, unless specified
  defval('g',fralmanac('GravAcc'));
  G=fralmanac('GravCst');
  defval('xver',0)
  if xver==1
    disp(sprintf('E= %5.3g; v= %5.3f',E,v))
  end

  % Wavenumbers are in radians per meter
  k=2*pi./lbd;

  % Turcotte and Schubert (3-72)
  % Flexural Rigidity [Pa*m^3 or N*m]
  D=(E*Te.^3)/(12*(1-v^2)); 

  % Create grid on which to calculate coherence-square
  if length(Te)>=1 && length(f2)==1 && length(r)==1
    [LL,DD]=meshgrid(lbd,D);
    FF2=f2; RR=r;
  elseif length(f2)>=1 && length(Te)==1 && length(r)==1
    [LL,FF2]=meshgrid(lbd,f2);
    DD=D; RR=r;
  elseif length(r)>=1 && length(Te)==1 && length(f2)==1
    [LL,RR]=meshgrid(lbd,r);
    DD=D; FF2=f2;
  else
    error('Only one vector allowed')
  end
  KK4=(2*pi./LL).^4;
  KK=(2*pi./LL);

  % Forsyth Eqs. (3) and (6)
  xai=1+DD.*KK4/drho/g;
  phi=1+DD.*KK4/rc/g;
  beta=rc./xai/drho;

  % Should you ever need to AVERAGE you would need to average the
  % elements in the numerator and the denominator separately  
  % Olhede and Simons, Eq. (64)
  Ctop=(xai+FF2.*rc^2./drho^2.*phi-RR.*sqrt(FF2)*rc./drho.*[phi.*xai+1]).^2;
  Cbot1=xai.^2+FF2.*rc^2./drho^2-2*RR.*sqrt(FF2)*rc/drho.*xai;
  Cbot2=1+FF2.*rc^2./drho^2.*phi.^2-2*RR.*sqrt(FF2)*rc/drho.*phi;
  % See Forsyth Eq. (25)
  G2b=Ctop./Cbot1./Cbot2;

  if xver==1
    % Simply put another way, Eq. (65)
    Ctop=(xai.*drho^2+FF2.*rc^2.*phi-RR.*sqrt(FF2)*rc.*drho.*[phi.*xai+1]).^2;
    Cbot1=xai.^2.*drho^2+FF2.*rc^2-2*RR.*sqrt(FF2)*rc*drho.*xai;
    Cbot2=drho^2+FF2.*rc^2.*phi.^2-2*RR.*sqrt(FF2)*rc*drho.*phi;
    difer(G2b-Ctop./Cbot1./Cbot2,6,[],NaN);
  end

  % Check the transition wavelength in TRANSL

  % Bouguer admittance in accord with Forsyth's Eqs (11)-(12)
  Zb=-2*pi*G*rc*exp(-KK*T).*...
     (1./xai+phi.*FF2.*beta.^2-RR.*sqrt(FF2)*rc./drho.*(phi.*xai+1)./xai.^2)./...
		   (1+FF2.*beta.^2-2*RR*sqrt(FF2)*rc./drho./xai);

  % Convert to mgal/m
  Zb=Zb/1e-5;

  % Free air admittance also in mgal/m
  Zf=Zb+2*pi*G*rc/1e-5;

  l=lbd;

  % Output
  varns={G2b,k,l,Zb,Zf};
  varargout=varns(1:nargout);
elseif strcmp(Te,'demo')
  % Illustrates the functions FORSYTH, MCKENZIE, ADMITTANCE, AND TRANSL
  % For EQUAL loading at both interfaces 
  % For UNCORRELATED loading only 
  % Elastic thickness [km]
  Te=[20 80];
  % Density contrasts in kg/m^3
  DEL=[2670,630];
  % Depth to the second loading interface [km]
  T=40;
  % Calculates coherence and admittance a la Forsyth
  [G2bF,k,l,ZbF,ZfF]=forsyth(Te*1e3,[],1,0,DEL(1),DEL(2),T*1e3);
  % Figures out the transitional wavelength
  [k12,l12]=transl(1,Te,DEL(1),DEL(2));
  % Calculates coherence and admittance a la McKenzie
  [l,ZbTe1,G2bTe1,ZfTe1]=mckenzie([0 DEL(1) DEL(1)+DEL(2)],[0 T],[1 1],Te(1));
  [l,ZbTe2,G2bTe2,ZfTe2]=mckenzie([0 DEL(1) DEL(1)+DEL(2)],[0 T],[1 1],Te(2));
  % Calculates Bouguer admittance a la old-fashioned version
  [Q1,kQ]=admittance(Te(1)*1e3,1,0,T*1e3,DEL(1),DEL(2));
  [Q2,kQ]=admittance(Te(2)*1e3,1,0,T*1e3,DEL(1),DEL(2));
  
  % Make the plot
  clf
  [ah,ha,H]=krijetem(subnum(3,1));
  axes(ah(1))
  pF=semilogx(k*1000,G2bF,'Color','g','LineW',2);
  hold on
  pMTe1=semilogx(2*pi./l*1000,G2bTe1,'b');
  pMTe2=semilogx(2*pi./l*1000,G2bTe2,'r');
  set([pMTe1 pMTe2],'MarkerS',4); 
  grid on; openup(gca,6); 
  xl(1)=xlabel('Wavenumber (rad/km)');
  yl(1)=ylabel('Bouguer Coherence \gamma^2');
  pl(1)=plot(k12(1),0.5,'x');
  pl(2)=plot(k12(2),0.5,'o');
  hold off

  axes(ah(2))
  pZf=semilogx(k*1000,ZfF*1000); hold on
  pZfTe1=semilogx(2*pi./l*1000,ZfTe1*1000,'b');
  pZfTe2=semilogx(2*pi./l*1000,ZfTe2*1000,'r');
  grid on; openup(gca,6); hold off
  xl(3)=xlabel('Wavenumber (rad/km)');
  yl(2)=ylabel('Free-air Admittance Z_f (mgal/km)');

  axes(ah(3))
  pZb=semilogx(k*1000,ZbF*1000); hold on
  pZbTe1=semilogx(2*pi./l*1000,ZbTe1*1000,'b');
  pZbTe2=semilogx(2*pi./l*1000,ZbTe2*1000,'r');
  pQ1=semilogx(kQ*1000,Q1*1000,'g');
  pQ2=semilogx(kQ*1000,Q2*1000,'y');
  grid on; openup(gca,6); hold off
  xl(4)=xlabel('Wavenumber (rad/km)');
  yl(3)=ylabel('Bouguer Admittance Z_b (mgal/km)');
  
  % Cosmetics
  set(ah,'xlim',[3e-3 1e-1])
  set(ah(1),'ytick',[0:0.25:1])
  xx(1)=xtraxis1d(ah(1)); 
  xl(2)=xlabel('Wavelength (km)');
  xx(2)=xtraxis1d(ah(2)); 
  xx(3)=xtraxis1d(ah(3)); 
  longticks([ah xx],2)
  fig2print(gcf,'tall')
  figdisp
end

function varargout=simulos(th0,params,xver)
% [Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb]=SIMULOS(th0,params,xver)
%
% Simulates data under the UNCORRELATED two-layer Forsyth model as used
% by Olhede & Simons with the primary spectra S11 etc being the
% equilibrium-loading ones. So we are not simulating initial loads and then
% performing the flexural elastic compensation; rather we are simulating
% the effects under such an elastic compensation model.
%
% S11 is the power spectral density of the EQUILIBRIUM load at interface 1.
%
% INPUT:
%
% th0      The true parameter vector with elements:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=s2   The first Matern parameter, aka sigma^2 
%          th(4)=nu   The second Matern parameter 
%          th(5)=rho  The third Matern parameter 
% params   A structure with constants that are (assumed to be) known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
%          z2    the positive depth to the second interface [m]
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs 0 Don't blur likelihood using the Fejer window
%                N Blur likelihood using the Fejer window [default: N=2]
%          kiso  [immaterial whether this is here or not]
%
% OUTPUT:
%
% Hx       Real matrix of spatial-domain observations [m], see Hk
% Gx       Real vector with Bouguer gravity anomaly data [mgal]
% th0      The true parameter vector pertaining to this data
% params   The structure with the known knowns, see above
% k        Wavenumber(s) suitable for the data sets returned [rad/m]
% Hk       A complex matrix of Fourier-domain observations, namely
%          final surface and subsurface topography [m]
% Gk       Complex vector with Fourier-domain Bouguer anomaly [m/s^2]
% Sb       The spectral matrix that you've used in this process
% Lb       The Cholesky decomposition of the spectral matrix which you
%          might use to evaluate the fit later on, in which case you
%          don't need any of the previous output
%
% EXAMPLE:
%
% simulos('demo1')
% simulos('demo2')
% simulos('demo3')
%
% SEE ALSO:
%
% MLEOS, LOADING, SIMULROS
%
% Last modified by fjsimons-at-alum.mit.edu, 04/30/2014

% Check how it behaves when NOT a power of two! FFT should still be exact
% so wouldn't matter. Implement the windowing!

% This corresponding to random guesses for the parameters in MLEOS
aguess=[1e24 0.8 0.0025 2 2e4];
%aguess=[7e22 0.4 0.0025 2 2e4];

% Here is the true parameter vector and the only variable that gets used 
defval('th0',aguess);

% Here is the extra verification parameters
defval('xver',1)

if ~isstr(th0)
  % Get or supply the needed additional parameters
  fields={'DEL','g','z2','dydx','NyNx','blurs'};
  defstruct('params',fields,...
	    {[2670 630],9.81,35000,[20 20]*1e3,[128 128]/2,2});

  if xver==1
    % Dump to screen
    osdisp(th0,params)
  end
  
  % Extract the variables explicitly from this structure
  for ind=1:length(fields)
    eval(sprintf('%s=params.(fields{ind});',fields{ind}))
  end
  
  % First make the wavenumbers, given the data size and the data length
  [k,dci,dcn]=knums(params);

  % This should make sense as the spacing in wavenumber domain
  dkxdky=2*pi./NyNx./dydx;

  % Now construct the whole-spectral matrix
  Z1=randgpn(k,dci,dcn);
  Z2=randgpn(k,dci,dcn);
  disp(sprintf('Z1: mean %+6.3f ; stdev %6.3f',...
	       mean(Z1(:)),std(Z1(:))))
  disp(sprintf('Z2: mean %+6.3f ; stdev %6.3f',...
	       mean(Z2(:)),std(Z2(:))))

  switch blurs
   case {0,1}
    disp('SIMULOS without BLURRING')
    % Now make the spectral-spectral portion of the spectral matrix
    S11=maternos(k,th0);
    % The Cholesky decomposition of the lithospheric-spectral matrix
    [~,~,L,T]=Tos(k,th0,params);
    % Roll in the sqrt of the factored portion
    Lb=repmat(sqrt(S11),1,3).*L;
    % The spectral matrix in case you care but you don't
    Sb=[S11.*T(:,1) S11.*T(:,2) S11.*T(:,3)];
   otherwise    
    % If I stay on the same k-grid, I'm really not doing any convolution at
    % all as you can see quickly. So by "suitably discretizing" the
    % convolutional operator we mean performing it on a highly densified
    % grid of which the target grid may be a subset. Doing this on the
    % same grid would be the "inverse crime" of not changing the grid at
    % all. Run Fk for this case to see it then would be a delta function
    
    % Blurs IS the refinement parameter; make new wavenumber grid
    disp(sprintf('SIMULROS with BLURRING factor %i',blurs))
    k2=knums(params,1);

    % Now make the spectral-spectral portion of the spectral matrix
    S11=maternos(k2,th0);
    % The lithospheric-spectral matrix on this second grid
    [~,~,~,T]=Tos(k2,th0,params); 
    % Which we multiply by the spectral-spectral portion
    S=[S11.*T(:,1) S11.*T(:,2) S11.*T(:,3)];
    
    % Now do the blurring and subsampling to original grid
    Sb=bluros(S,params,1);
    
    % imagesc(decibel(Fejk))
    % Old, isotropic, unblurred
    % imagesc(decibel(reshape(S(:,1),NyNx)))
    % New, blurred, crossy on new grid
    % imagesc(decibel(reshape(Sb(:,1),blurs*NyNx)))
    % New, blurred, crossy, subsampled on old grid
    % imagesc(decibel(conv2(Fejk,reshape(S(:,1),blurs*NyNx),'same')))
    
    % And then we do the Cholesky decomposition of that, explicitly
    Lb=[sqrt(Sb(:,1)) Sb(:,2)./sqrt(Sb(:,1)) ...
	sqrt((Sb(:,1).*Sb(:,3)-Sb(:,2).^2)./Sb(:,1))];

    % Could do CHOLCHECK here
    if xver==1
      cholcheck(Lb,Sb,6,1)
      cholcheck(Lb,Sb,6,2)
    end
        
    % Should make sure that this is real! 
    Lb=realize(Lb);
  end
  % Blurred or unblurred, go on

  % And put it all together, unwrapped over k and over x
  Hk=[Lb(:,1).*Z1(:) [Lb(:,2).*Z1(:)+Lb(:,3).*Z2(:)]];
  
  % Without the L we should probably get the equilibrium topographies
  % So for Airy icebergs we should just have a simple scaling?
  if th0(1)==0 && th0(2)==0
    % This should be 1
    difer(Lb(:,1)-1,[],[],NaN)
    % This should have some relation to airyratio
    airyratio=params.DEL(1)/params.DEL(2);
    difer(Lb(:,2)+airyratio,[],[],NaN)
    % This should be 0
    difer(Lb(:,3),[],[],NaN)
  end
    
  % Now create the gravity observation also
  G=fralmanac('GravCst');
  % With this sign convention the depth is positive
  Gk=2*pi*G*DEL(2)*exp(-k(:).*z2).*Hk(:,2);

  % And go to the space domain - unitary transform
  Hx(:,1)=tospace(Hk(:,1),NyNx);
  Hx(:,2)=tospace(Hk(:,2),NyNx);
  Gx     =tospace(Gk     ,NyNx);

  % Return the output if requested
  varns={Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb};
  varargout=varns(1:nargout);
elseif strcmp(th0,'demo1')
  svnm=th0;
  [Hx,Gx,th0,p,k,Hk,Gk]=simulos;
  clf
  kelicol
  [ah,ha]=krijetem(subnum(2,2));
  [tl(1),cb(1),xc(1),xa(1)]=plotit(ah(1),Hx(:,1)/1000,size(k),...
		       'final surface','topography (%s)','km');
  [tl(2),cb(2),xc(2),xa(2)]=plotit(ah(2),Hx(:,2)/1000,size(k),sprintf(...
      'final subsurface z_2 = %i km',round(p.z2/1000)),'topography (%s)','km');
  [tl(3),cb(3),xc(3),xa(3)]=plotit(ah(3),Gx*1e5,size(k),sprintf(...
      'gravity anomaly with %s%s = %i kg m^{-3}',...
      '\Delta','\rho',p.DEL(2)),'Bouguer anomaly (%s)','mgal');

  % Maybe we can calculate the coherence here also?
  %   NW=3;
  %   [FX,FY,SX,SY,SXY,COH2]=...
  %     mtm3(reshape(Hx(:,1),NyNx),reshape(Gx,NyNx),NW,max(round(2*NW)-2,1));
  
  NyNx=p.NyNx;
  dydx=p.dydx;

  % Cosmetics
  they=linspace(1,NyNx(1),5);
  thex=linspace(1,NyNx(2),5);
  spunkm=(NyNx-1).*dydx/1000;
  set(ah,'ylim',they([1 end])+[-1 1]/2,...
	 'xlim',thex([1 end])+[-1 1]/2,...
	 'ytick',they,...
	 'xtick',thex,...
	 'ytickl',-spunkm(1)/2+(they-1)*dydx(1)/1000,...
	 'xtickl',-spunkm(2)/2+(thex-1)*dydx(2)/1000)
  longticks([ah cb])
  nolabels(ah(1:2),1)
  nolabels(ha(3:4),2)
  
  % Plot the parameters here
  axes(ah(4))
  axis([-1 1 -1 1])
  nolabels(ah(4)); noticks(ah(4)); box on; axis image
  
  xof=-0.75;
  tx(1)=text(xof, 0.75,sprintf('%s = %12.3g','D',       th0(1)));
  tx(2)=text(xof, 0.50,sprintf('%s = %12.3g','f^2',     th0(2)));

  tx(3)=text(xof, 0.00,sprintf('%s = %12.3g','\sigma^2',th0(end-2)));
  tx(4)=text(xof,-0.25,sprintf('%s = %12.3g','\nu',     th0(end-1)));
  tx(5)=text(xof,-0.50,sprintf('%s = %12.3g','\rho',    th0(end)));
  moveh(ah(4),-getpos(cb(2),3))
  
  fig2print(gcf,'portrait')
  figdisp([],svnm,[],1,'pdf')
elseif strcmp(th0,'demo2')
  % Try to get close to the example we had in 2000
  NyNx=[100 100];
  D=2.4722e+23;
  f2=0.5;
  s2=0.001;
  nu=0.5;
  rho=30000;
  z2=15000;
  % Perform the simulation
  [Hx,Gx,th0,k,Hk,Gk,params]=simulos;
elseif strcmp(th0,'demo3')
  params.blurs=2;
  params.NyNx=[128 128];
  [Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb]=simulos([],params); 
  Lk=Lkos(k,th0,params,Hk); 
  imagesc(decibel(v2s(Lk)))
  axis image
end

% Plotting routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tls,cbs,xcb,xa]=plotit(aha,dats,nm,stronk,strink,unid)
axes(aha)
imagesc(reshape(dats,nm)); 
axis image
limc=halverange(dats,95,NaN);
% Later, be more sophisticated than this
if limc(1)>-1
  limc=round(100*limc)/100;
else
  limc=round(limc);
end

%shrink(aha,1.2,1.2)
set(aha,'clim',limc)
tls=title(stronk);
xa=xlabel(sprintf('mean %+6.3f ; stdev %6.3f %s',...
		  mean(dats(:)),std(dats(:)),unid));
cbs=colorbar('ver');
axes(cbs)
xcb=ylabel(sprintf(strink,unid));
set(cbs,'ylim',limc)

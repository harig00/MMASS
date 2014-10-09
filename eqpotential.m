function varargout=eqpotential(M,deplonlat,Lup,Sstrain,catflag)
% [phiE1,phiE1r,Sstrain]=EQPOTENTIAL(M,[dep lon lat],Lup,Sstrain,catflag)
%
% Computes the first-order Eulerian gravity potential perturbation due to
% an earthquake, in a spherical Earth, at the surface, as real spherical
% harmonic expansion coefficients, as 'lmcosi', suitable for expansion
% with PLM2XYZ. Note that if you use PLM2POT you have to pull out the
% GM/r as this function will be putting it in. Also note that DELTAJL
% will pull out an extra factor of sqrt(2*pi) as they use unit normalized
% harmonics and we don't, here.
%
% INPUT:
%
% M            Moment tensor [Mrr Mtt Mpp Mrt Mrp Mtp] [default: vertical strike slip]
%              Units are in the customary dyne cm, OR:
%              An earthquake CMT solution string
% dep/lon/lat  Earthquake depth/longitude/latitude (km, degrees) [default: 15,180,0]
%              Note that depth is with respect to the radius, and there
%              may be an ocean layer on top of the model!, OR:
%              A path string where a CMT-solution file will be found
% Lup          Truncation degree of the ROTATED coefficients [default: none]
% Sstrain      Optional structure array containing nn,el,ww,a,b,c,d,Pre,Urfs
%              the combined outputs from GETSPHEROIDAL and SMODESTRAIN
% catflag      0 Modes files made by YANNOS aka MINOSY [default]
%              1 Modes files made by MINOS_BRAN and formatted by
%                EIGCON from  the packages MINEOS-1.0.0 or MINEOS-1.0.2; 
%                the former is a small catalog, the latter a big
%                one... but remember the Earth models are slightly different!
%              2 Modes files in the format read by LOADMODES
%
% OUTPUT:
%
% phiE1        Degree/order/cosine/sine coeffficients for earthquake at North Pole
% phiE1r       Degree/order/cosine/sine coeffficients at true earthquake location
%              Units are properly those of potential, J/kg or m^2/s^2
% Sstrain      Structure array with nn,el,ww,a,b,c,d,Pre,Urfs used here
%
% EXAMPLES:
%
% eqpotential('demo1') % Recreate something like Dahlen and Tromp (1998) Table 5.1 
% eqpotential('demo2') % Kuril Islands 2006/11/15 (M200611151114A)
% eqpotential('demo3') % Kuril Islands 2007/11/15
% eqpotential('demo4',0) % 2010 Chilean Earthquake with YANNOS
% eqpotential('demo4',1) % 2010 Chilean Earthquake with MINEOS-1.0.0
%
% SEE ALSO: GETSPHEROIDAL, PLM2XYZ, PLM2ROT, SMODESTRAIN, MODESUM
%
% Last modified by efwelch-at-princeton.edu, 07/21/2010
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% Supply defaults
defval('M',[1 1 -2 0 0 0]/sqrt(6))
defval('Sstrain',NaN)
defval('catflag',0)
defval('xver',0)

if ~isstr(M) || isempty(strmatch('demo',M))
  if ~isstr(M)
    % M was a moment tensor thus we supply more numbers
    defval('deplonlat',[15 180 0])
  else
    % M was an CMT event identifier thus we supply a path
    defval('deplonlat',fullfile(getenv('IFILES'),'CMT','jan76_feb10.ndk'));
    % Now read the file CMTSOLUTION-M
    M=cmtsol(M,deplonlat);
    % Assign to the right variables
    deplonlat=[M.Dep M.Lon M.Lat];
    name=M.EventName;
    M=M.M.*10^M.Exp;
  end

  if xver==1
    % Calculate my moment, see SMOMENT
    % Mw=2/3/log(10)*...
    %        log(1/sqrt(2)*sqrt(M(1)^2+M(2)^2+M(3)^2+...
    % 			  2*(M(4)^2)+2*(M(5)^2)+2*(M(6)^2)))-10.7;  
    Mw=smoment(M);
    disp(sprintf('Moment magnitude Mw = %3.1f',Mw))
  end
  
  % Convert from [dyne cm] to Joules [J]
  M=M*1e-7;
  % Get depth, longitude and latitude; convert depth to m
  dep=deplonlat(1)*1000;
  % Convert from geographic to geocentric coordinates, see GEOCENTRICRTHETA
  lon=deplonlat(2);
  lat=deplonlat(3);
  
  % Display what's being done
  if xver==1
    disp(sprintf('Moment tensor is:'))
    disp(sprintf('| %9.2g %9.2g %9.2g |',M(1),M(4),M(5)))
    disp(sprintf('| %9.2g %9.2g %9.2g |',M(4),M(2),M(6)))
    disp(sprintf('| %9.2g %9.2g %9.2g |',M(5),M(6),M(3)))
  end
  
  % If not provided, get the outputs from a defaulted GETSPHEROIDAL
  % If adding arguments, making Sstrain NaN will generate it, too
  if nargin<4 || [~isstruct(Sstrain) && isnan(Sstrain)]
    clear Sstrain
    [rad,nn,el,ww,U,V,P,dUdr,dVdr,dPdr]=getspheroidal([],[],catflag);
    % Select the excitation coefficient corresponding to the source depth
    rs=rad(end)-dep;
    Sstrain=smodestrain(rs,[],[],rad,nn,el,ww,U,V,P,dUdr,dVdr);
    clear rad nn el ww U V P dUdr dVdr dPdr
  end

  % Make sure the right flag was set
  difer(min(size(Sstrain.Pre))-1,[],[],NaN)

  % The order-independent prefactors for the mode summation
  stuff=Sstrain.Pre(:)'.*Sstrain.ww(:)'.^(-2).*((2*Sstrain.el(:)'+1)/4/pi);

  % Note that the accuracy needs to be fixed by delaying the
  % normalization in GETSPHEROIDAL etc, see the latest unincorporated
  % changes... the Sstrain.Pre are the prefactors at the hypocenter
  
  % Perform the normal mode summation over the branches
  a=modesum(Sstrain.a.*stuff,Sstrain.nn,Sstrain.el);
  b=modesum(Sstrain.b.*stuff,Sstrain.nn,Sstrain.el);
  c=modesum(Sstrain.c.*stuff,Sstrain.nn,Sstrain.el);
  d=modesum(Sstrain.d.*stuff,Sstrain.nn,Sstrain.el);

  if xver==1
    disp(sprintf('Maximum degree in catalog is %i',max(Sstrain.el)))
  end
  
  % Check these all have equal length
  alle=[length(a) length(b) length(c) length(d)];
  difer(sum(diff(alle)),[],[],NaN); difer(sum(alle)-4-4*max(Sstrain.el),[],[],NaN)
  alle=alle(1);

  % Now, for a typical expansion, need to know where these shall go
  [dems,dels,mz,phiE1]=addmon(max(Sstrain.el));

  % Normalization. The above expressions are for straight Legendre
  % cosine and sine expansions using unnormalized associated Legendre
  % fuctions a la (10.53). However, our PLM2XYZ expands into 4\pi
  % normalized real spherical harmonics without the Condon-Shortley
  % phase... in other words, must prep the coefficients for this.
  % The factor, however, is different depending on the order 0<=m<=2.
  els=[0:dels(end)]';
  k=sqrt(els.*(els+1));
  fax0=1./sqrt(2*els+1);
  fax1=fax0/sqrt(2).*k;
  fax2=fax1.*sqrt(k.^2-2);

  % Multiply the moment tensor at this point - orders are separate
  % Stick in these elements - remember that the sequences start from l=m
  % DT (10.54)-(10.59)
  phiE1(dems==0,3)=fax0       .*(a*M(1)+c*[M(2)+M(3)]);
  phiE1(dems==1,3)=fax1(2:end).*(d(2:end)*M(4));
  phiE1(dems==1,4)=fax1(2:end).*(d(2:end)*M(5));
  phiE1(dems==2,3)=fax2(3:end).*(b(3:end)*[M(2)-M(3)]);
  phiE1(dems==2,4)=fax2(3:end).*(2*b(3:end)*M(6));

  % Do a basic check of the multipolarity
  if xver==1
    if difer(phiE1(dems==0,3:4)) ; disp(sprintf('Monopole')) ; end
    if difer(phiE1(dems==1,3:4)) ; disp(sprintf('Dipole')); end
    if difer(phiE1(dems==2,3:4)) ; disp(sprintf('Quadrupole')); end
  end
  
  % Save time in what's next by truncating the series
  defval('Lup',dels(end))

  if Lup~=dels(end)
    if Lup>0
      disp(sprintf('Rotated coefficients truncated to L = %i',Lup))
    end
  end

  if Lup>0
    % Then rotate the man to the required location
    % Change is to now say 180-lon
    phiE1r=plm2rot(phiE1(1:addmup(Lup),:),0,90-lat,180-lon);
    % Don't forget to take take a quick look at the spectrum
  else
    phiE1r=phiE1;
  end
  
  % Provide output
  varn={phiE1,phiE1r,Sstrain};
  varargout=varn(1:nargout);
elseif strcmp(M,'demo1')
  % Compare with Dahlen and Tromp (1998) Table 5.1 
  clear M
  % Need to make this a bit better of course
  % Need to make a polar plot... evaluate at pole, without rotating
  % Now check at least the pattern
  M{1,1}=[ 1  1  1  0 0 0 ]/sqrt(3); 
  M{2,1}=[ 0  0  0  0 0 1 ]/-sqrt(2); 
  M{3,1}=[ 0  0  0  1 0 0 ]/sqrt(2); 
  M{4,1}=[ 1 -1  0  0 0 0 ]/sqrt(2); 
  M{5,1}=[ 1  1 -2  0 0 0 ]/sqrt(6);  
  M{6,1}=[-2  1  1  0 0 0 ]/sqrt(6); 
  M{1,2}=[ 1  1  1  0 0 0 ]/-sqrt(3); 
  M{2,2}=[ 0  1 -1  0 0 0 ]/sqrt(2); 
  M{3,2}=[ 0  0  0  0 1 0 ]/sqrt(2); 
  M{4,2}=[ 1  0 -1  0 0 0 ]/sqrt(2); 
  M{5,2}=[ 1 -2  1  0 0 0 ]/sqrt(6); 
  M{6,2}=[-2  1  1  0 0 0 ]/-sqrt(6); 

  % Produce fictitious earthquake at the equator
  defval('deplonlat',[16 180 0])

  % Begin figure
  clf
  [ah,ha]=krijetem(subnum(6,2));
  % Perhaps should focus on the maximum degree of GRACE as anything else
  % will not be detected anyways.
  Lup=359; degres=0.5;
  Lup=120; degres=0.5;
  Lup=24; degres=0.5;
  Lup=60; degres=0.5;
  % Actually, save time, do NOT rotate by setting Lupr to 0
  Lupr=Lup;
  % Make the polar grid if we aren't rotating
  if Lupr==0
    % The size could relate to degres... if you knew r already
    [phi,theta]=polargrid(2*round(180/degres)+1,4*round(180/degres)+1);
  end
  
  plottype=1;
  switch plottype
   case 1
    jays=1:2;
   case 2
    jays=1;
  end

  % Do the calculation
  for i=1:6
    for j=jays
      if i==1 & j==1
	[p,pr,Sstrain]=eqpotential(M{i,j},deplonlat,Lupr);  
      else
	[p,pr]=eqpotential(M{i,j},deplonlat,Lupr,Sstrain);
      end
      % Display
      ix=(i-1)*2+j;
      axes(ah(ix))
      switch plottype
       case 1
	% Expand the result - up to Lup
	[r,lon,lat]=plm2xyz(pr(1:addmup(Lup),:),degres);
	if Lupr==0
	  % Rather than rotating, interpolate on a polar grid
	  % Remember we're in essence flipping the dimensions
	  % Can make it even faster by interpolating only the center points
	  % rather than the whole set and then truncating the axes.
	  ri=interp2(lon/180*pi,(90-lat)/180*pi,r,phi,theta);
	  % Plot on the map
	  imagesc(ri); axis image;
	  caxis([-max(abs(ri(:))) max(abs(ri(:)))]/2)
	  axis(([size(ri,2) size(ri,2) size(ri,1) size(ri,1)]-1)/2+...
	       [-1 1 -1 1]*length(ri)/60)
	  % Is it worth checking the contour lines? Use CONTOURC?
	  hold on ; [c,h]=contour(ri,[0 0]); set(h,'LineC','k')
	  % Is it worth comparing the plot with PLOTPLM option 5?
	else
	  % Plot rotated function on the map
	  imagesc(r); axis image;
	  caxis([-max(abs(r(:))) max(abs(r(:)))]/2)
	  axis(([size(r,2) size(r,2) size(r,1) size(r,1)]-1)/2+...
	       [-1 1 -1 1]*length(r)/60)
	  hold on ; [c,h]=contour(r,[0 0]); set(h,'LineC','k')
	end
	cb(ix)=colorbar('vertical');
       case 2
	plot(3:p(end,1),decibel(indeks(plm2spec(p),'4:end')))
	xlim([0 p(end,1)])
	% yim([-100 0])
      end
    end
  end
  % Cosmetics
  fig2print(gcf,'tall')
  figdisp
  if plottype==1
    kelicol
    nolabels(ah)
    t=supertit(ah(1:2),sprintf('Depth = %i km ; maximum degree L = %i',...
			       round((Sstrain.re-Sstrain.rs)/1000),Lup));
  else
    longticks(ah,2)
    nolabels(ha(1:5),1)
    xlabel('degree l')
    ylabel('power (dB)')
    delete(ha(7:end))
  end
elseif strcmp(M,'demo2')
  % KURIL ISLANDS 2006/11/15  M200611151114A
  M=[1.740  -0.555   -1.180   1.640   2.580  -0.771]*1e28;
  % This the PDE location
  deplonlat=[38.9 153.29 46.57];
  % This the CENTROID location
  deplonlat=[13.5 154.33 46.71];
  % Make the prediction
  [phiE1,phiE1r]=eqpotential(M,deplonlat);
  % Now you can also call it directly with the CMT codename
  % M='M200611151114A';
  % [phiE12,phiE1r2]=eqpotential(M);
  % difer(phiE1-phiE12)
  % difer(phiE1r-phiE1r2)
  data=plotplm(phiE1r,[],[],4,[]);
  TH=5;
  [lon,lat]=caploc([deplonlat(2) deplonlat(3)],5);
  hold on
  plot(lon,lat,'k'); axis image
  axis([102 217 15 73])
  % Now convert the potential perturbation to Slepian coefficients
  L=60;
  % To make this easier on GLMALPHAPTO, round the angles for reuse
  % The "round"s didn't use to be there
  [falpha,V,N,MTAP,C]=plm2slep(phiE1r,TH,L,...
			       round(deplonlat(2)),...
			       90-round(deplonlat(3)),0);
  % Now got to do something with this don't I
  [i,j]=sort(abs(falpha),'descend');
  % plot(falpha(j(1:30)),'o-')
elseif strcmp(M,'demo3')
  % KURIL ISLANDS 2007/11/15
  M=[-1.380   1.330  0.051  -0.274  -0.785   0.798];
  deplonlat=[10.0  154.52 46.24];
  % Make the prediction
  [phiE1,phiE1r]=eqpotential(M,deplonlat);
  plotplm(phiE1r,[],[],4,[])
elseif strcmp(M,'demo4')
  keep=M;
  M='C201002270634A'; % 2010 Chilean earthquake
  %M='M122604A'; % 2010 Sumatran earthquake
  Mstr=M;
  
  % Which mode catalog? The SECOND argument in the demo
  defval('deplonlat',[]);
  catflag=deplonlat;
  defval('catflag',1);

  % Load the event
  M=cmtsol(M);
  deplonlat=[M.Dep M.Lon M.Lat];
  M=M.M*10^M.Exp;
  
  % Axis range for the plot (with respect to the epicenter)
  lora=10; lara=10;

  % See if you've already done it, and save it
  fdir=fullfile(getenv('IFILES'),'EQPOTENTIAL');
  fname=sprintf('%s/%s_%i.mat',fdir,'chilefagrav',catflag);
  %fname=sprintf('%s/%s_%i.mat',fdir,'sumatrafagrav',catflag);
  if exist(fname,'file')==2 & 1==3
    load(fname)
    % Due to an earlier version the saved file had the degree 2
    % coefficient of the predicted coseismic potential anomaly taken
    % out. I suppose I should rerun and resave. It's PLM2POT that changed. 
    disp(sprintf('%s loaded',fname))
  else
    % Make the prediction and save the modes
    [phiE1,phiE1r,S]=eqpotential(M,deplonlat,[],[],catflag);
    % Do remember that PLM2POT assumes that the potential needs to be
    % multiplied by GM/r
    defval('GM',fralmanac('GM_EGM96','Earth'))
    defval('r',fralmanac('a_EGM96','Earth'))
    % Note this is always in here, don't confuse when using PLM2XYZ
    phiE1r(:,3:4)=phiE1r(:,3:4)./GM*r;
    % Change potential to free-air gravity anomaly, expand to Lup
    Lup=62;
    % Reference to nothing since it's all a perturbation
    d=plm2xyz(plm2pot(phiE1r(1:addmup(Lup),:),[],[],[],2,'nothing'));
    % Convert from m/s^2 to microgal
    d=d*1e8;
    save(fname,'d','M','deplonlat','catflag','S','phiE1r')
  end  

  % Make the plot
  clf
  fig2print(gcf,'portrait')
  ah=gca;
  
  opsjun=2;
  if opsjun==1
    cax=round(halverange(d));
    [h,cax]=imagefnan([0 90],[360 -90],d,[],cax);
    plotcont
    xli=deplonlat(2)+[-lora lora];
    yli=deplonlat(3)+[-lara lara];
    axis([xli+[xli<0]*360 yli]); 
    [cb,cbs]=addcb([],cax,cax,'kelicol',range(cax)/10);
    movev(cb,-0.035)
    set(cbs,'string',sprintf('free-air anomaly (%sgal)','\mu'))
    longticks([cb],2)
    longticks([ah cb],2)
    deggies(ah)
    figdisp([],sprintf('%s_%i',keep,catflag))
  elseif opsjun==2
    % Full-resolution of EGM2008
    L=1500;
    degres=180/sqrt(L*(L+1));

    % Redo the expansion for a more complete suite, see above
    % By using plm2pot we simply take out the (2,0) term instead of
    % really referring it to WGS84, but the effect should be small
    Lup=60;
    % Reference to nothing since it's all a perturbation
    r=plm2xyz(plm2pot(plmfilt(phiE1r,Lup),[],[],[],2,'nothing'),degres);
    % Plot it
    clf
    ah=gca;
    [r,c,ph]=plotplm(setnans(r));

    kelicol
    tt=kelicol;
    tt(81,:)=[1 1 1];
    tt(82,:)=[1 1 1];
    colormap(tt)

    % A scale as in PLM2DIFF for the monthlies
    caxis([-1 1]*1e2*1e-9*1.25)

    cb=colorbar('hor');
    axes(cb)
    xlabel(sprintf(...
        'free-air anomaly to L = %i predicted due to %s [m/s^2]',Lup,Mstr))
    set(cb,'xaxisl','t')

    shrink(cb,2,2)
    axes(cb)
    longticks(cb,2)
    movev(cb,-.05*2)

    fig2print(gcf,'portrait')
    figdisp

    % Make it bigger to get a good bounding box
    set(ah,'camerav',6)
    movev(cb,-.05*2)
    movev([ah cb],.1)

    moveh([ah cb],-.015)
    %print -djpeg -r300 /home/fjsimons/EPS/alumnitea3_300
    %print -djpeg -r600 /home/fjsimons/EPS/alumnitea3_600
    print -djpeg -r300 /home/fjsimons/EPS/alumnitea4_300
    print -djpeg -r600 /home/fjsimons/EPS/alumnitea4_600
  end
end

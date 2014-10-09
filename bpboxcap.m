function varargout=bpboxcap(TH,Lmax,meth,evens,sord,xver)
% [Bl,dels]=BPBOXCAP(TH,Lmax,meth,evens,sord,xver)
%
% Returns the spherical harmonic power spectrum of cylindrical 
% boxcar polar caps: 1 in the region and 0 elsewhere.
% Dahlen and Simons (2008) eqs (48-50)
%
% INPUT:
%
% TH       Colatitudinal radius of the cap, in degrees, may be vector
%          Colatitudinal halfwidth of the belt, degrees, may be vector
% Lmax     Maximum degree of this non-bandlimited function
% meth     1 Analytically, using unnormalized Legendre polynomials [default]
%          2 Numerically, integrating normalized Legendre polynomials 
% evens    0 For all degrees [default]
%          1 Only return the even degrees [logical for equatorial symmetry]
% sord     1 Single cap of diameter 2TH [default]
%          2 Double cap left by subtracting belt of width 2TH
%          3 Equatorial belt of width 2TH
% xver     1 Provides extra verification
%          0 Doesn't [default]
%
% OUTPUT:
%
% Bl       The power spectrum 
% dels     The spherical harmonic degrees
%
% EXAMPLE:
%
% bpboxcap('demo1') % should return no messages and pass all checks
% bpboxcap('demo2') % compare with the results of GEOBOXCAP
%
% SEE ALSO:
%
% GEOBOXCAP, GAMMAP
%
% Last modified by fjsimons-at-alum.mit.edu, 02/03/2010

% For all the others, sum((2*dels+1).*Bl)=4*pi does NOT hold since they
% are not normalized that way.

if ~isstr(TH)
  defval('TH',100)
  defval('Lmax',100)
  defval('meth',1)
  defval('sord',1)
  defval('xver',0)
  
  t=cputime;
  
  % If double boxcar cap, always symmetric, always faster, reinsert at the end
  if sord==2 || sord==3
    waseven=evens;
    defval('evens',1);
  else
    defval('evens',0)
  end
  
  % Returns the spherical harmonic degrees
  dels=0:evens+1:Lmax;
  dels1=0:Lmax+1;
  switch meth 
   case 1
    if sord==1
      % Get the unnormalized spherical harmonics
      [pel,mu,norms]=plm(dels1,0,cos(TH*pi/180),xver);
    else
      [pel,mu,norms]=plm(dels1,0,sin(TH*pi/180),xver);
    end
    % Get the power spectrum by application of the analytical formula
    Bl=pi*([repmat(1,1,length(TH)) ; pel(1+evens:evens+1:end-2,:)]...
	   -[pel(2:evens+1:end,:)]).^2 ...
       ./repmat((2*dels'+1).^2,1,length(TH));
    if sord==2 || sord==3
      Bl=Bl.*[1+(-1).^repmat(dels',1,length(TH))].^2;
    end  
    if sord==3 % Need extra term for the zeroth order
	       % Here you cannot take out the prefactors
	       Bl(1,:)=pi*(2-2*[repmat(1,1,length(TH))-pel(2,:)]).^2;
    end
   case 2
    switch sord
     case 1
      intv=[cos(TH(:)*pi/180) repmat(1,length(TH),1)];
     case 2
      intv=[cos(pi/2-TH(:)*pi/180) repmat(1,length(TH),1)];
     case 3
      intv=[cos(pi/2+TH(:)*pi/180) cos(pi/2-TH(:)*pi/180)];
    end
    % Overdo the GL weights for all but the last one
    [wGL,xGL]=gausslegendrecof(Lmax,[],intv);
    % Initialize this matrix as a potentially 3D object
    xel=repmat(NaN,[length(dels) size(xGL)]);
    % Do Gauss-Legendre on a whole matrix at the same time
    % Writing (:,:,:) explicitly helps when there is only one degree
    xel(:,:,:)=xlm(dels,0,acos(xGL),xver);
    % Initialize the results matrix
    Bl=repmat(NaN,length(dels),length(TH));
    for index=1:length(TH)
      % Get the power spectrum by numerical integration
      Bl(:,index)=...
	  [(1+[sord==2])*2*pi*wGL(:,index)'*xel(:,:,index)'].^2./(2*dels+1);
    end
  end
  
  if (sord==2 || sord==3) 
    if evens==1 && waseven==0
      Bll=zeros(length(dels1)-1,length(TH));
      dels=dels1(1:end-1);
      Bll(1:2:end,:)=Bl;
      Bl=Bll;
    elseif evens==0 && waseven==0
      % Strictly speaking, for method 2 this is only necessary for sord=2
      Bl(2:2:end,:)=0;
    end
  end
  
  goods=sprintf(' BPBOXCAP A: %6.3f',sqrt(Bl(1,:)*4*pi));
  % Check the B0 term which should equal the area^2 divided by 4pi in the
  % unit-normalized basis, where the Y00 term equals 1/sqrt(4*pi)
  A=4*pi*spharea(TH,sord);
  difer(Bl(1,:)-A.^2/4/pi,[],[],goods)
  
  % Prepare output
  varns={Bl,dels};
  varargout=varns(1:nargout);

elseif strcmp(TH,'demo1')
  [b1,l1]=bpboxcap([20 30 40],50,1,1,1);
  [b2,l2]=bpboxcap([20 30 40],50,2,1,1); difer(b1-b2)
  [b1,l1]=bpboxcap([20 30 40],50,1,0,1);
  [b2,l2]=bpboxcap([20 30 40],50,2,0,1); difer(b1-b2)
  [b1,l1]=bpboxcap([20 30 40],50,1,1,2);
  [b2,l2]=bpboxcap([20 30 40],50,2,1,2); difer(b1-b2)
  [b1,l1]=bpboxcap([20 30 40],50,1,0,2);
  [b2,l2]=bpboxcap([20 30 40],50,2,0,2); difer(b1-b2)
  [b1,l1]=bpboxcap([20 30 40],50,1,0,3);
  [b2,l2]=bpboxcap([20 30 40],50,2,0,3); difer(b1-b2)
  [b1,l1]=bpboxcap([20 30 40],50,1,1,3);
  [b2,l2]=bpboxcap([20 30 40],50,2,1,3); difer(b1-b2)
  b=bpboxcap(180,ceil(rand*200),1);
  difer(2*sum(b)-b(1)-4*pi)
elseif strcmp(TH,'demo2')
  % Make a circle somewhere on the equator
  N=500;
  TH=ceil(rand*60);
  lon1=180+ceil((rand-0.5)*90);
  lat1=0+ceil((rand-0.5)*90);
  [lon2,lat2]=caploc([lon1 lat1],TH,N,1);
  % Make a picture also
  clf; [ah,ha]=krijetem(subnum(2,1));
  axes(ah(1)); plot(lon2,lat2,'k'); hold on; plot(lon1,lat1,'ko')
  title(sprintf('TH = %i ; %s_0 = %i ; %s_0 = %i',...
		TH,'\phi',lon1,'\theta',90-lat1))
  axis([0 360 -90 90]); grid on
  set(gca,'xtick',unique([0 lon1 360]),'ytick',unique([-90 lat1 90]))
  xver=0; meth=1+1*round(rand); meth=1;
  Lmax=90;
  % Calculate the power spectral density of a boxcar of the same size
  [Bl,dels]=bpboxcap(TH,Lmax,meth,0,1,xver);
  % Calculate the same for this particular boxcar using another method
  % For this here, smaller is better
  degres=0.5;
  [Bl2,dels2]=geoboxcap(Lmax,[lon2 lat2],N,degres,1);
  axes(ah(2))
  semilogy(dels,Bl,'b+-')
  hold on
  semilogy(dels2,Bl2,'rv')
  title(sprintf('rmse = %6.3e',rms(Bl-Bl2)))
end


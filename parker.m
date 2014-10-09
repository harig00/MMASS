function tercor=parker(c11,cmn,topo,ninterf,addm,drho)
% tercor=PARKER(c11,cmn,topo,ninterf,addm,drho)
% 
% Computes terrain/Bouguer corrections using Parker's method.
% Bathymetry input is negative. The resulting correction must be ADDED to
% the free-air anomaly. Everything in standard units (SI)!
% So 'tercor' is in m/s2. To convert to milliGal, divide by 1e-5.
% Take any nans' out, as in
% tdint(isnan(tdint))=mean(tdint(~isnan(tdint)));
% See Parker, GJRAS 1972 (31) 447-455
% Compare this to adding 2*pi*drho* the water depth to "replace"
% the water with mantle density.
%
% INPUT:
%
% c11,cmn        Coordinates of top left and bottom right pixel [degrees]       
% topo           Topography/bathymetry matrix [m],
% ninterf        Number of interfaces
% addm           Depth to the interface(s) in [m>0]
% drho           Density contrast(s) [kg/m3]
%
% tercor=parker(c11,cmn,topo,2,[0 6000],[1700 600]);
%
% Last modified by fjsimons-at-alum.mit.edu, October 23rd, 2003
% Contributions from Mark Behn, Ban-Yuan Kuo and Jian Lin

%---Define Input Parameters---
npower=5;            % Power of Taylor Series Expansion
ifold=1;             % Folding flag (ifold=1 : folding, ifold=2 : no folding)

%---Load Datafile---
[ny,nx]=size(topo);
lenx=abs(cmn(1)-c11(1))*fralmanac('DegDis');
leny=abs(c11(2)-cmn(2))*fralmanac('DegDis');
dx=lenx/(nx-1);
dy=leny/(ny-1);

%---Normalize data to new reference level-----
zmax=max(topo(:));
zmin=min(topo(:));
slev=mean(topo(:));
ddepth=slev-topo; %$ How much deeper than mean depth

% Gravitational constant in SI
grav=2*pi*drho*fralmanac('GravCst'); 

if ifold==1 
  ddepth=flipflop(ddepth);
  lenx=2*lenx;
  leny=2*leny;
end

% Wave numbers rad/m
K=knum2(size(ddepth),[leny lenx]);
% Because fft is not fftshifted here not either
kwn=fftshift(K);

% Parker (1972) Eq. 4
csum=zeros(size(ddepth));
for ip=1:npower
  fprintf('%5s  %12.4f\n','Power =',ip)
  data=ddepth.^ip;
  data=fft2(data);
  data=data.*(kwn.^(ip-1))/factorial(ip);
  csum=csum+data;
end

[k,l]=find(kwn==0);
csum(k,l)=0;

fprintf('%10s\n','Upward Continuation')
tercor=zeros(ny,nx);

for iref=1:ninterf
  % Depth over which to be upward continued, starting from mean elevation
  zlev=abs(slev)+addm(iref);
  % Upward continuation: the field must be attenuated
  data=csum*grav(iref).*exp(-zlev*kwn); 
  
  fprintf('%10s\n','Inverse Transformation')
  data=ifft2(data);
  if mean(abs(imag(data(:))))>1e-14
    warning(sprintf('Transformed data is not real by %s',...
		    mean(abs(imag(data(:))))))
  end
  ddepth=real(data(1:ny,1:nx));
  tercor=tercor+ddepth;
end


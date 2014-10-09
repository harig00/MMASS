function [th,k12Mth,kMth,g2M,mlazi,Mlazi,xyml,xyMl,a,b]=halfmazi(data,minxmaxx,varargin)
% [th,k12Mth,kMth,g2M,mlazi,Mlazi,xyml,xyMl,a,b]=halfmazi(data,[minx maxx])
% [th,k12Mth,kMth,g2M,mlazi,Mlazi,xyml,xyMl,a,b]=halfmazi(data,[minx maxx],[cenx ceny],thin)
%
% thin          Return this particular theta(thin)
%
% OUTPUT
%
% th            Azimuthal angles, in degrees
% k12Mth        Wavenumber of half-maximal coherence as a function of theta
% kMth          Wavenumber of the full maximum as a function of theta
% g2M           Maximum value of the coherence on the line of this azimuth
% mlazi         The theta of mazimum k12Mth (WEAKEST direction)
% Mlaz          The theta of minimum k12Mth (STRONGEST direction)
% xyml          Plots a line through the coherence at the STRONGEST direction
% xyMl          Plots a line through the coherence at the WEAKEST direction
%               (yes, the capital letters are defined backwards here)
% a             Wavenumbers and  
% b             Interpolated coherence at the requested theta(thin)
%
% For a two-dimensional matrix 'data' calculates, in function 
% of a line of constant azimuth going to the center of it,
% the location of the maximum value along this line,
% as well as the location of the half maximum. For symmetric square
% matrices, if needed give the symmetry center explicitly.
% Put in the smallest non-zero wavenumber and the largest one (assuming
% the plot is symmetric with the same wavelength ranges in the x and y
% direction.)
% 
% See also HALFMISO, the isotropic equivalent.
%
% CHECK IF THE CENTERS ARE RIGHT BY HALTING THE PROGRAM AND
% PLOTTING THINGS UP.
%
% Check the precision parameter for the wavelength dropoff.

% Last modified by fjsimons-at-mit.edu, November 30th, 2005
% Re: the precision parameter.

precirat=1/40;

[ry,rx]=size(data);
hrx=ceil(rx/2);
hry=ceil(ry/2);

if hrx~=hry
  error('Should input square data set')
end

[minx,maxx]=deal(minxmaxx(:,1),minxmaxx(:,2));

resol=180;
th=linspace(0,180,resol);
if nargin>=3
  senter=varargin{1};
else
  senter=[hrx hry];
end
if nargin>=4
  thin=varargin{2};
end

for index=1:length(th)
  % Work with the minimum number of independent measurements
  % for a given azimuth; assume the center chosen is the symmetry center!
  kset=min(size(data));
  [x,y]=linazim(th(index),senter,hrx,[],kset);
  % Get the coherence data, interpolated on the line through the center
  % at the indexed azimuth the
  polated=interp2(data,x,y,'linear');
  % Since it is symmetric, only get the second half of this line
  polpos=polated(senter(1):end);
  % Length of the original section through the coherence
  orlen=1:kset-senter(1)+1;
  % Requested length, at even more points
  intlen=linspace(1,kset-senter(1)+1,10*kset);  
  warning off
  % Interpolate that coherence curve through many more points
  polat2=interp1(orlen,polpos,intlen,'cubic');
  warning on
  % Find the maximum of this curve
  [g2M(index),kMthIJ]=max(polat2);
  % Find where coherence first drops to preci*half maximum.
  preci=g2M(index)*precirat;
  howclose=abs(polat2-g2M(index)/2);
  k12MthIJ=min(find(howclose<preci));
  if ~isempty(k12MthIJ)
    % Find the first wavenumber that satisfies this drop-off
    k12Mth(index)=minx+(intlen(k12MthIJ)-1)/(intlen(end)-1)*(maxx-minx);
    kMth(index)=minx+(intlen(kMthIJ)-1)/(intlen(end)-1)*(maxx-minx);
  else
    k12Mth(index)=NaN;
    kMth(index)=NaN;
  end
  [aj,callpro]=star69;
  if strcmp(callpro,'australiadir2') | strcmp(callpro,'australiadir3') | ...
	 strcmp(callpro,'australiadir40') | strcmp(callpro,'jonathan') 
    if index==thin
      a=linspace(minx,maxx,length(intlen));
      b=polat2;
    end
    if index==1
      clf
      subplot(211)
      % Gets fed flipud(data); plot up again
      % imagef([-maxx  maxx],[maxx -maxx],flipud(data))
      imagesc(data)
      axis image; hold on
      longticks      
      axis xy
      nolabels      
    end
    if index==1 | ~mod(index,20)  
      subplot(211)
      plot(x,y,'w-')
      subplot(212)
      %      plot(linspace(minx,maxx,length(orlen)),polpos,'r+'); 
      plot(linspace(minx,maxx,length(intlen)),polat2,'b-'); 
      xlim([minx maxx])
      hold on ; set(gca,'xscale','log')
      if ~isnan(k12Mth(index))
	sp1=plot(k12Mth(index),polat2(k12MthIJ),'o');
	%      sp2=plot(kMth(index),g2M(index),'o'); 
	xlo=xlim;
	%      plot(xlim,[g2M(index)/2 g2M(index)/2],'--'); 
	%      plot(xlim,[g2M(index)/2 g2M(index)/2]+preci,'--','Color',[1 1 1]*0.7); 
	%      plot(xlim,[g2M(index)/2 g2M(index)/2]-preci,'--','Color',[1 1 1]*0.7); 
	set([sp1],'MarkerF','w','MarkerE','k','MarkerS',8)
      end    
    end    
  end
end

% Maximum wavenumber, minimum wavelength, thus weakest theta
Mlazi=indeks(th(k12Mth==min(k12Mth)),1);
mlazi=indeks(th(k12Mth==max(k12Mth)),1);
[xml,yml]=linazim(Mlazi,senter,hrx,[],min(size(data)));
[xMl,yMl]=linazim(mlazi,senter,hrx,[],min(size(data)));
xyml=[xml(:) yml(:)];
xyMl=[xMl(:) yMl(:)];

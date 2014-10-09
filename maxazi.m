function [azim,pol,azimax,nav,x,y]=maxazi(data,varargin)
% [azim,pol,azimax,nav,x,y]=MAXAZI(data)
% [azim,pol,azimax,nav,x,y]=MAXAZI(data,azim)
% [azim,pol,azimax,nav,x,y]=MAXAZI(data,azim,[cenx ceny])
% [azim,pol,azimax,nav,x,y]=MAXAZI(data,azim,[cenx ceny],[minx maxx])
%
% OUTPUT:
%
% azim        Azimuthal directions
% pol         Average coherence along these directions
% azimax      Direction of the maximum average coherence
% nav         Number of coherence samples in this average
% x,y         The indices through the coherence that contains the maximum
% 
% May specify center of plot if it is not in the geometric center of the
% indices. For FFT products hence symmetric about the origin!
% If a minimum and maximum wavenumber is specified, will compare
% the azimuthal curves with LOGLINEARLY spaced samples. In other words,
% the mean won't be biased as badly towards the many short-wavelength
% coherence values.
%
% Calculates mean 'pol' and maximum 'azimax' of mean azimuthal line integral
% for a data set. Or just calculates the integral over the given azimuth
% Mean calculated on the basis of 'nav' samples, from 0 to 180.
% Suitable for kx, ky wavelength representation: limit to kx^2+ky^2 rings.
% Also returns the x and y indices used to define the maximum azimuth.
% 
% See also ISAV, RADAVG, MINAZI, HALFMAZI
%
% Note that for images, may need to do flipud before calculating.
%
% EXAMPLE:
%
% maxazi('demo');
%
% Last modified by fjsimons-at-mit.edu, November 30th, 2005

if ~isstr(data)
  [ry,rx]=size(data);
  hrx=ceil(rx/2);
  hry=ceil(ry/2);
  
  if hrx~=hry
    error('Should input square data set')
  end

  resol=180;
  if nargin==1 | isempty(varargin{1})
    azim=linspace(0,180,resol);
  else
    azim=varargin{1};
  end
  if nargin>=3
    senter=varargin{2};
  else
    senter=[hrx hry];
  end
  if nargin==4
    minx=varargin{3}(1);
    maxx=varargin{3}(2);
    disp('Loglinear mean comparison!')
  else
    disp('Linear mean comparison!')
  end

  for index=1:length(azim)
    % Work with the minimum number of independent measurements
    % for a given azimuth; assume the center chosen is the symmetry center!
    % 'senter' and min(size(data)) dictate if indeed senter(1):end is a good
    % choice (see e.g. halfmazi)
    
    [x,y]=linazim(azim(index),senter,hrx,[],min(size(data)));
    polated=interp2(data,x,y,'linear');
    if nargin==4
      logpol=interp1(linspace(minx,maxx,length(polated(senter(1):end))),...
		     polated(senter(1):end),logspace(log10(minx),log10(maxx),...
						     length(polated((senter(1):end)))));
      % Compare
      % plot(linspace(minx,maxx,length(polated(senter(1):end))),...
      %	polated(senter(1):end),'+-'); hold on
      % plot(logspace(log10(minx),log10(maxx),length(polated(senter(1):end))),...
      %  	logpol,'ro-')
      pol(index)=nanmean(logpol);  
      nav(index)=length(logpol);
    else
      pol(index)=nanmean(polated(senter(1):end));  
      nav(index)=length(polated(senter(1):end));
    end
  end

  [i,k]=max(pol);  
  azimax=azim(k);

  [x,y]=linazim(azimax,senter,hrx,[],min(size(data)));
elseif strmatch(data,'demo')
  [X,Y]=meshgrid(linspace(-10,10,200));
  G1=gauss2(X,Y,3,6);
  [azim,pol]=maxazi(G1);
  subplot(221); imagesc(G1); axis image
  subplot(222) ; plot(azim,pol); xlim([0 180]); grid
  [Xrot,Yrot]=rotdom(linspace(-10,10,200),linspace(-10,10,200),pi/6);
  G2=gauss2(Xrot,Yrot,3,6);
  [azim,pol]=maxazi(G2);
  subplot(223); imagesc(G2); axis image
  subplot(224) ; plot(azim,pol); xlim([0 180]); grid
end

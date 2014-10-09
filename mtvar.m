function varargout=mtvar(Sl,l,L,TH,sord)
% [varSl,Sl,l]=MTVAR(Sl,l,L,TH,sord)
%
% After calculating a multitaper power spectral density, calculates the
% appropriate error bars based on it, according to Dahlen & Simons
% (2008). Not for whole-sphere data, for which the simple formula in 
% Dahlen & Simons eq. (47) applies, see also RB VII p. 32.
%
% INPUT:
%
% Sl     The power spectral density whose variance we are calculating
% l      Degrees at which the Sl is being quoted [lengths must match]
% L      Bandwidth of the tapers considered [if 1, whole-sphere]
% TH     String identifying a region, e.g. 'australia' (sord=0), OR
%        Colatitudinal radius of the cap, in degrees <=180 (sord=1), OR
%        Colatitudinal halfwidth of the cut, degrees <=90 (sord=2)
%        FJS TO FIX TO ALLOW XY CURVE ALSO
% sord   1 Single cap of diameter 2TH [default]
%        2 Double cap left by subtracting belt of width 2TH
%        0 Automatically supplied in case TH is a region string 
%
% OUTPUT:
%
% varSl  The variance of the power spectral density
% Sl     The power spectral density again, should you require it
% l      The spherical harmonic degrees again, if you want them
%
% EXAMPLES:
%
% mtvar('demo1') % Comparison with MVARRATIOS
% mtvar('demo2') % Comparison as TH grows very large, which should tend
%                  to one over 2L+1 or so, as per our paper
% mtvar('demo3') % Comparison between analytic and prescribed circle
%
% SEE ALSO: MVARRATIOS, MULTIVAR, GAMMAP, GAMMAB, MCOVARIANCES
%
% Last modified by fjsimons-at-alum.mit.edu, 10/04/2010

% Default is the effect on a white spectrum
defval('l',[0:100])
defval('Sl',repmat(1,1,length(l)))

if ~isstr(Sl)
  if length(Sl)~=length(l)
    % Any combination of [Sl l] is allowed, no requirement on the degrees
    % to be adjacent or complete, calculation only needs a single pair 
    error('Number of degrees must match length of the power spectrum')
  end
  
  % Go the extra mile for verification, or not
  defval('xver',0)
  defval('TH','australia')
  defval('L',8)

  % Figure out what sort of region we're doing here
  if isstr(TH)
    % Override whatever may have been given
    sord=0;
  else
    defval('sord',1);
  end

  if L==1
    disp('Whole-sphere estimate')
    varSl=2./[2*l+1].*Sl.^2;
  else
    % Highest degree of a zero-j database, else, switch to approximation
    zjmax=500;
    
    % Highest degree of a six-j database, else, switch to scaled full-sphere
    % or method 3 for the GAMMAP calculation
    sjmax=40;
    
    % The maximum degree used from the zero-j database under this construction
    Lmax=min(max(l),zjmax);
    
    % Get all the ZEROJ coefficients at the same time
    [allW,C0,S0,Leff]=zeroj(repmat(0:2:2*L,1,Lmax+1),...
			    gamini(0:Lmax,L+1),gamini(0:Lmax,L+1));
    
    if 2*L<=sjmax
      % Always only get the evens since we're studying l=l, the variance
      [Gp,p,K]=gammap(L,TH,sord,1,1);
      disp(sprintf('Shannon number  %5.2f',K))
      A=spharea(TH,sord);
      disp(sprintf('Fractional area %5.2f%s',100*A,'%'))
    else
      error('Fix using the methods outlined in MVARRATIOS')
    end
    
    % Initialize the varSl variance ratio
    varSl=repmat(NaN,1,length(l));
    
    % Step through the required degrees
    for ixl=1:length(l)
      if l(ixl)<=Lmax
	% Stick with the one-blow precalculated ones if it's below the bandwidth
	W=allW((L+1)*l(ixl)+1:(L+1)*(l(ixl)+1));
      else
	% Calculate on the fly rather than from the loaded database
	W=wigner0j(2*L,l(ixl),l(ixl),xver);
	W=W(1:2:end);
	if xver==1
	  % Compare to the asymptotic representation if you want 
	  md=sum(abs(W-(-1)^l(ixl)*...
		     [plm(0:2:2*L,0,0)./sqrt((2*l(ixl)+1))]')/(L+1));
	end
      end
      % Only select the evens since we're doing VARIANCE at equal l=l'
      % This is DS eq. (165)
      varSl(ixl)=Sl(ixl)^2./(2*pi)*[(2*p+1).*Gp]*[W.^2]';
    end
  end
  
  % Prepare output
  varns={varSl,Sl,l};
  varargout=varns(1:nargout);

elseif strcmp(Sl,'demo1')
  disp('Comparison with the results of MVARRATIOS')
  L=round(rand*20);
  TH=10+round(rand*30);
  sord=1+round(rand);
  disp(sprintf('\nChecking result for L = %i ; TH = %i ; sord = %i\n',...
	       L,TH,sord))
  l=1:100;
  % Calculate the ratio the old way (same routine though, just checking)
  [mt2ws,l,mt2wsinf]=mvarratios(L,TH,sord,l);
  % Calculate the absolute variance the new way
  Sl=ones(size(l));
  [varSl,Sl,l]=mtvar(Sl,l,L,TH,sord);
  % Make a plot of this variance
  clf ; errorbar(l,Sl,sqrt(varSl))
  title(sprintf('\nL = %i ; TH = %i ; sord = %i\n',L,TH,sord))
  axis tight
  % We should be getting the result that this here is zero
  difer(varSl/2.*(2*l+1)-mt2ws)
elseif strcmp(Sl,'demo2')
  keyboard
end


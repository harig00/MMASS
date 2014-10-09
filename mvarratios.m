function [mt2ws,l,mt2wsinf]=mvarratios(L,TH,sord,l,xver)
% [mt2ws,l,mt2wsinf]=MVARRATIOS(L,TH,sord,l,xver)
%
% Calculates the ratio of the variance of the eigenvalue-weighted
% multitaper estimates to the untapered whole-sphere estimation variance 
% which we can derive from first-principle statistics...
% Dahlen and Simons (2008) eq. (178)
%
% INPUT:
%
% L            Bandwidth of the multitaper windows, always a scalar
% TH           Colatitudinal radius of the cap, in degrees <=180, may be vector
%              Colatitudinal halfwidth of the cut, degrees <=90, may be vector
% sord         1 Single cap of diameter 2TH [default], always a scalar
%              2 Double cap left by subtracting belt of width 2TH, always scalar
% l            Degrees at which the calculation is to be carried out 
%
% OUTPUT:
%
% mt2ws        Ratio of the estimation variance of the eigenvalue-weighted
%              multitaper estimate to the untapered whole-sphere estimate
%              This is a matrix of dimensions length(TH) x length(l)
% l            The degrees at which this variance ratio was computed [0->50]
% mt2wsinf     Asymptotic approximation to the above for degrees >> bandwidth
%
% SEE ALSO: PERIODOVAR, MULTIVAR, MTVAR
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/29/2007

% SHOULD USE THIS FUNCTION IN DSVARRATIO3

% Supply default values
defval('L',20)
defval('TH',30)
defval('sord',1)
defval('l',0:50)
defval('xver',1)

% Highest degree of a zero-j database, after that, switch to approximation
zjmax=500;
% Highest degree of a six-j database, after that, switch to scaled
% full-sphere result - actually, with the improvements, best switch to
% method 3 if you can do this
sjmax=40;
% The maximum degree
Lmax=min(max(l),zjmax);

% Get all the ZEROJ coefficients at the same time
[allW,C0,S0,Leff]=zeroj(repmat(0:2:2*L,1,Lmax+1),...
			gamini(0:Lmax,L+1),gamini(0:Lmax,L+1));

% NEED TO BUILD IN WHAT HAPPENS WHEN A IS EQUAL TO ZERO

% Calculate the area under the taper
A=4*pi*spharea(TH,sord);
if sord==3
  error('This choice is not allowed')
end

if 2*L<=sjmax
  % Calculate the matrix that goes into this - this takes most of the time
  % Always only get the evens since we're studying l=l, the variance
  [Gp,p,K]=gammap(L,TH,sord,1,1);
  disp(sprintf('We are getting a Shannon number of %5.2f',K))
  needscaling=0;
elseif 2*L>sjmax
  % Cannot do the full works for lack of the 6j coefficients
  % But instead simply scale the A=4pi multitaper result
  % Initialize some important arrays
  bigS=gamini([0:L],(L+1))';
  bigSp=repmat([0:L]',(L+1),1);
  cG=repmat(NaN,2*L+1,1);
  % Calculate the curly Gamma instead of the regular gamma
  for e=0:2*L
    cG(e+1,1)=sum([2*bigS'+1].*[2*bigSp'+1].*...
		  zeroj(bigS,e,bigSp,Leff,[],C0,S0).^2);
  end
  % Turns out that this "curly Gamma" is what regular Gamma goes to when A
  % goes to the whole sphere (A=4pi) --- compare only the even degrees
  % Under the single cap thing: single cap of entire globe
  Gp=cG(1:2:end)'*4*pi/(L+1)^4;
  p=0:2:2*L;
  % After this, no sense in verifying since the verification wants to
  % know the truth, which you know you can't know
  xver=0;
  % This part of the loop was verified to give the same results even for
  % small L as intended
  needscaling=1;
end

% Initialize the mt2ws variance ratio
mt2ws=repmat(NaN,length(TH),length(l));

% Better get all of the wigner0j symbols at once here
for ixl=1:length(l)
  if l(ixl)<=Lmax
    % Stick with the one-blow precalculated ones if it's below the bandwidth
    W=allW((L+1)*l(ixl)+1:(L+1)*(l(ixl)+1));
  else
    % Calculate on the fly rather than from the loaded database
    W=wigner0j(2*L,l(ixl),l(ixl),xver);
    W=W(1:2:end);
    % Compare to the asymptotic representation if you want 
    md=sum(abs(W-(-1)^l(ixl)*[plm(0:2:2*L,0,0)./sqrt((2*l(ixl)+1))]')/(L+1));
    %   disp(sprintf(...
    %	'With the asymptotic form you''d have errors of %5.3e but luckily' ...
    %	' you didn''t'],md))
  end
  % Only select the evens since we're doing VARIANCE at equal l=l'
  mt2ws(:,ixl)=(2*l(ixl)+1)/(4*pi)*[repmat(2*p+1,length(TH),1).*Gp]*[W.^2]';
end
% Now use scaling if you need to to rectify
sexp=-0.88;
if needscaling==1
  mt2ws=mt2ws*(A/4/pi)^sexp;
  disp(sprintf('MVARRATIOS Using scaling (A/4%s)^%5.3f','\pi',sexp))
end

% FOR SMALL L VERIFY THAT THIS APPROACH GIVES REASONABLE RESULTS BOTH FOR
% THE SINGLE AND THE DOUBLE CAP NOTE THAT K MUST BE VERY BIG
% mvarratios(20,10,2,0:50) % is the double cap with cut 10
% Now come in here and hand-edit to see that it works.

% What should the large-l limit be? Use the Legendre functions instead
% Initialize the large-l limit of this
mt2wsinf=[1/(4*pi)*[repmat(2*p+1,length(TH),1).*Gp]*...
	  plm(0:2:2*L,0,0).^2]';
%disp(sprintf(['MVARRATIOS: Error of l=inf approximation at'...
%	      ' l/L=%3i is %4.1e%s'],l(end)/L,...
%	      abs(mt2ws(end)-mt2wsinf*((A/4/pi)^(sexp*(needscaling==1))))/...
%	      mt2ws(end)*100,'%'))

if xver==1 && l(1)==0 && ~[sord==1 & TH==0]
  % What should the variance ratio be at l=0 for each of these TH's?
  % We have three ways to calculate this
  mt2ws01=repmat(NaN,1,length(TH));
  mt2ws02=repmat(NaN,1,length(TH));
  mt2ws03=repmat(NaN,1,length(TH));
  % Initialize some important arrays
  bigS=gamini([0:L],(L+1))';
  bigSp=repmat([0:L]',(L+1),1);
  cG=repmat(NaN,2*L+1,1);
  % Calculate the "curly Gamma"
  for e=0:2*L
    cG(e+1,1)=sum([2*bigS'+1].*[2*bigSp'+1].*...
		  zeroj(bigS,e,bigSp,Leff,[],C0,S0).^2);
  end
  % Turns out that this curly Gamma is what regular Gamma goes to when A
  % goes to the whole sphere (A=4pi) --- compare only the even degrees
  % Under the single cap thing: single cap of entire globe
  difer(cG(1:2:end)'*4*pi/(L+1)^4-gammap(L,180,1,1,1),[],[],...
	'Check for A=4pi from single cap passed')
  % Under the double cap thing: subtract belt of nothing
  difer(cG(1:2:end)'*4*pi/(L+1)^4-gammap(L,0,2,1,1),[],[],...
	'Check for A=4pi from double cap passed')
  for inx=1:length(TH)
    % Use the power spectrum of the single/double cap
    [BeL(inx,:),dels]=bpboxcap(TH(inx),2*L,[],0,sord);
    % Compute this alternative
    mt2ws01(inx)=1/4/pi/K(inx)^2*[(2*[0:2*L]+1).*BeL(inx,:)]*cG;
  end
  % Compare just for good measure
  difer(mt2ws01'-mt2ws(:,1),[],[],...
	'MVARRATIOS: First check for zero-l passed')
  % Another check for the l=0 ratio
  for inx=1:length(TH)
    [G,V,EL,EM,KN]=glmalpha(90*(sord==2)+(-1)^(sord+1)*TH(inx),L,sord);
    mt2ws02(inx)=sum(V.^2)/sum(V)^2;
  end
  difer(mt2ws02'-mt2ws(:,1),[],[],...
	'MVARRATIOS: Secnd check for zero-l passed')
  % There is a third check which we should put in there as well
  % See also Wieczorek and Simons 2005!
  if sord==1
    for inx=1:length(TH)
      % Get the square localization matrix for axisymmetric domains
      D=dlmlmp(TH(inx),L,L,xver);
      % Have since rewritten this for the double cap, as well
      mt2ws03(inx)=sum(abs(D(:).^2))/K(inx).^2;
    end
    difer(mt2ws03'-mt2ws(:,1),[],[],...
	  'MVARRATIOS: Third check for zero-l passed')
  end
end


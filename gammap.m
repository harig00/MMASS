function [Gp,p,K]=gammap(L,TH,sord,evens,meth,isk,xver)
% [Gp,p,K]=GAMMAP(L,TH,sord,evens,meth,isk,xver)
%
% Returns the funny GAMMA function from Dahlen & Simons (2008) eq. (167)
% For axisymmetric (order-degenerate functions), this affords a
% huge computational savings for methods 2 and 3 which are largely
% unworkable otherwise.
% Remember, the tapers are normalized to 4pi and not 1.
%
% INPUT:
%
% L         Taper bandwidth that is required in this case, always a scalar
% TH        Colatitudinal radius of the cap, in degrees <=180, may be vector
%           Colatitudinal halfwidth of the cut, degrees <=90, may be vector
%           ... OR: the Shannon number... see option 'isk'
%           ... OR: a string identifying a region, e.g. 'australia'
% sord      1 Single cap of diameter 2TH [default]
%           2 Double cap left by subtracting belt of width 2TH
%           0 Automatically supplied in case TH is a region string
% evens     1 only the even p's [default, good for variance but not covariance]
%           0 all p's 
% meth      1 using Wigner 6j and boxcar spectra [fast, default]
%           2 using GL integration Gamma_p(alpha,beta) [slow]
%           3 method using Wigner 3j symbols Gamma_p(alpha,beta) [less slow]
%           i.e. they compute Gammap for all taper combinations and then
%           take the eigenvalue-weighted average... rather slow still
% isk       0 nothing unusual here [default]
%           1 the TH is actually a Shannon number K, and TH is calculated
%           3 allowable option when TH is a region string
% xver      1 Excessive verification and loop testing
%           0 None of that
% 
% OUTPUT:
%
% Gp        Gamma(p) as defined by Dahlen & Simons (2008), length(TH) x length(p)
% p         The degrees at which this is evaluated, containing evens or not
% K         The Shannon number(s)
%
% EXAMPLES:
%
% gammap('demo1') % Tests single-cap results and should return no errors
% gammap('demo2') % Tests double-cap results and should return no errors
%
% SEE ALSO:
%
% GAMMAB which employs the same algorithm but saves the (alpha,beta) parts
%
% NOTE:
%
% GLMALPHA and GALPHA have different conventions for TH - change them
% (t)here. Certain efficiency gains are still possible...
%
% difer(isnan(gammap(10,0,1))-isnan(gammap(10,90,2))) % No cap, all cut
% Check the last demo, put in the A=0 and A=4*pi cases explicitly
%
% Last modified by fjsimons-at-alum.mit.edu, 05/10/2011

% MUST MAKE SURE THAT THE ZERO-A RESULT IS A SPECIAL CASE HERE
% PERHAPS MUST MAKE SURE THAT THE FULL-SPHERE RESULT IS A SPECIAL CASE HERE

if ~isstr(L)
defval('L',10)
defval('TH',30)
defval('sord',1)
defval('evens',1)
defval('meth',1)
defval('isk',0)
defval('xver',0)

% The dimension of the basis
Lpot=(L+1)^2;

% Initialize arrays
if isstr(TH)
  % Only one string region
  lTH=1;
  % Reset any other values that might have been called for
  isk=3;
  sord=0;
else
  lTH=length(TH);
end
GpB=repmat(NaN,lTH,2*L+1);
% Select the evens only, perhaps
GpB=GpB(:,1:evens+1:end);
KB=repmat(NaN,lTH,1);

% For (multiple) integer THs or for region strings
ndto=0;
for ind=1:lTH
  switch isk
    case 0
     fnpl=sprintf(['%s/%i/GAMMAP-%i-%i-%i-%i-%i.mat'],...
		  fullfile(getenv('IFILES'),'GAMMAP','LTH'),...
		  L,TH(ind),L,sord,evens,meth);
     case 1
     fnpl=sprintf(['%s/GAMMAP-%i-%i-%i-%i-%i.mat'],...
		  fullfile(getenv('IFILES'),'GAMMAP','LK'),...
		  TH(ind),L,sord,evens,meth);
   case 3
     fnpl=sprintf(['%s/GAMMAP-%s-%i-%i-%i.mat'],...
		  fullfile(getenv('IFILES'),'GAMMAP','REGIONS'),...
		  TH,L,evens,meth);
   otherwise
    error('Specify valid option for isk')
  end

  if exist(fnpl,'file')==2 %& 1==3 % Force to recalculate, maybe
    % Got these already
    load(fnpl)
    disp(sprintf('Loading %s',fnpl))
    GpB(ind,:)=Gp;
    KB(ind,:)=K;
  else
    disp(sprintf('%s does not exist',fnpl))
    % Still need to do these here
    ndto=ndto+1;
    if ~isstr(TH)
      THdo(ndto)=TH(ind);
      lTHdo=length(THdo);
    else
      % For the regions, only one at a time, so far
      THdo=TH;
      lTHdo=1;
    end
  end
end

% Still need to do a bunch here; do them all at once with method 1
% or do a single additional one using method 2 or 3
if ndto>=1
  if ndto>1 && meth~=1
    error(sprintf('Method %i only for (additional) scalar TH',meth))
  end
  % Don't do these - the missing ones have NaN's in them
  GpBndo=GpB;
  KBndo=KB;
  % Number of (non-)zero integrals performed - if zero, could have saved
  % the time, clearly... although inside the loop we catch them early, too
  numnz=0;  numz=0;
  % necessary

  switch isk
   case 3
    % Calculate the area from the input region name TH
    A=4*pi*spharea(THdo);
    % Calculate the Shannon numbers at once
    KB=Lpot*A/(4*pi);
   case 0
    % Calculate the area from the input angle TH
    A=4*pi*spharea(THdo,sord);
    if sord==3
      error('This choice of sord is not allowed')
    end
    % Calculate the Shannon numbers at once
    KB=Lpot*A/(4*pi);
   case 1
    KB=THdo;
    % Calculate the area from the input Shannon number TH
    A=KB*4*pi/Lpot;

    switch sord % Watch out here for the definitions in GALPHA and GLMALPHA
      % And now give up the travesty and call THdo the angles
     case 1
      THdo=acos(1-A/2/pi)*180/pi;
     case 2
      THdo=asin(1-A/4/pi)*180/pi;
     otherwise
      error('This choice of sord is not allowed')
    end
    % Watch out for formally correct but meaningless results
    KB=KB(~imag(THdo) & THdo>=0);
    THdo=THdo(~imag(THdo) & THdo>=0);
    % If it doesn't make sense, the code returns NaN and only the sensible ones
    % will be saved. What if both empty? Do nothing
    if isempty(THdo) && isempty(KB)
      % Return whatever you got if you don't need to do any more
      Gp=GpB;
      K=KB;
      disp('Don''t bother with illegal A/K/L combinations')
      return
    end
  end
  
  % Re-initialize after possible deletion of certain imaginary THdo as before
  GpB=repmat(NaN,lTHdo,2*L+1);
  
  if ~isstr(THdo)
    % Get the power spectrum of the boxcar cap, force to return only
    % evens when there is equatorial symmetry
    [Be,ee]=bpboxcap(THdo,2*L,[],sord-1,sord);
    % If it's a DOUBLE cap, will need/get only even degrees e
  else
    % If it's a GEOGRAPHIC region
    [Be,ee]=geoboxcap(2*L,THdo);
  end
  
  % What are the combinations of 4 indices each of which varies between 0
  % and L? Let's see - make this and then apply the selection rules
  % Make a big matrix index vector and then apply the rules - 
  % or construct one in a for loop and apply the rules at every step if
  % this is more memory efficient
  % May use NDGRID for this! 
  bigS=        gamini([0:L],(L+1)^3)';
  bigSp=repmat(gamini([0:L],(L+1)^2)',(L+1)^1,1);
  bigU= repmat(gamini([0:L], L+1)',   (L+1)^2,1);
  bigUp=repmat(       [0:L]',         (L+1)^3,1);
  
  % Now we do a loop on the index p from zero to 2*L, maybe only on the evens
  % Now we do a loop on the index e from zero to 2*L, maybe only on the evens
  % The degrees at which this is evaluated
  % Only require even p's when the variance is calculated at a degree l,
  % not the covariance; this is for the output
  % Only require even e's when the double cap is calculated, for internal
  % use only
  
  if meth~=1
    % If it's a GEOGRAPHIC region, will need to rewrite this!
    % Get the coefficients for the tapers - you always need the
    % eigenvalues
    % Watch the conventions of GLMALPHA and GRUNBAUM2 have different conventions
    [GLMAB,V,ELG,EMG]=glmalpha(90*(sord==2)+(-1)^(sord+1)*THdo,L,sord,0);
    % Ahem, hold on, we are REDEFINING these functions HERE in this paper
    newdef='4*pi';
    disp(sprintf('GAMMAP Redefined normalization for tapers to %s',newdef))
    GLMAB=GLMAB*sqrt(eval(newdef));
    % Get some other things too
    [EM,EL,mz,blkm]=addmout(L);
    % Block sorted array for future reference 
    % FJS ELGB=EL(blkm);
    EMGB=EM(blkm);
  end

  if meth~=2
    % Must load EXACTLY the one with the specified bandwidth, no default
    % substitutions are to be accepted, if not, won't know which vector
    % you've got!!! So this is how to do this
    dbL=2*L;
    % Load the ZEROJ database just once, for the maximum at all required
    [~,C0,S0]=zeroj(0,0,0,dbL);
    keyboard
    % Load the THREEJ database just once, as well
    [~,C3,S3]=threej(0,0,0,0,0,0,dbL);
  end

  if meth==1
    % What if you don't have the ones between 30 and 40, as is the case?
    if dbL>30
      if dbL>40
	error(['Making a 6j database for degrees > 40 is a no-no: change' ...
	       ' meth to 3'])
      else
	% Use the precomputed L=40 data bases
	dbL=40;
      end
    end
    % Load the SIXJ database just once, for the maximum at all required
    [~,C6,S6]=sixj(0,0,0,0,0,0,dbL);
  end
  
  if meth==3
    % Don't forget to make them complex for the Wigner 3j formalism to work
    U=ummp(L);
    % Do it! Note this doesn't work for block-sorted arrays
    GLMAB=U*GLMAB;
    % The important thing is to do this only once, as well...
    % This is actually completely unnecessary, may as well skip it.
  end
    
  if xver==1
    disp('GAMMAP: verifying all we can, Sir')
  end
  
  disp(sprintf('GAMMAP Using METHOD %i',meth))

  % For every desired degree
  for p=0:evens+1:2*L
    disp(sprintf('Working on %i of %i',p,2*L))
    switch meth
     case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % This is the fancy 6-j method
      dsum=0;
      % Check something else here - for any p, Research Book VI page 5
      if xver==1
	dver=0;
	for e=0:2*L
	  % The next think is a general rule that must have ALL e's not
	  % just the even ones, it has nothing to do with single/double caps
	  disp('Excessive verification of some important things here')
	  % Note: the 2*L supplied must be exactly the same as the loaded
          % database! Don't go around changing this in any way.
	  % This is the blue hatched region in my notebook
	  dver=dver+...
	       (-1)^(p+e)*(2*e+1)*...
	       sixj(bigS,repmat(e,size(bigS)),bigSp,...
		    bigU,repmat(p,size(bigS)),bigUp,...
		    dbL,[],C6,S6)...
	       .*zeroj(bigS,e,bigSp,dbL,[],C0,S0).*...
	       zeroj(bigU,e,bigUp,dbL,[],C0,S0);
	end
	difer(dver-zeroj(bigS,p,bigUp,dbL,[],C0,S0)...
	      .*zeroj(bigU,p,bigSp,dbL,[],C0,S0))
      end
      
      % Note that the next loop definition only works for ROWS!!
      for e=ee(:)'
	% Check something else here - for p=0
	if xver==1 && p==0
	  % Gamma_0 needed when l is zero
	  % This is the green hatched region in my notebook
	  difer((L+1)^4-...
		sum((2*bigS+1).*(2*bigSp+1).*(2*bigU+1).*(2*bigUp+1).*...
		    [zeroj(bigS,p,bigUp,dbL,[],C0,S0).^2.*...
		     zeroj(bigSp,p,bigUp,dbL,[],C0,S0).^2]'))
	end
	
	% Selection rules - those of the 3j contain those of the 6j
	selx=triangle(bigS,p,bigUp) & triangle(bigSp,p,bigU) ...
	     & triangle(bigS,e,bigSp) & triangle(bigU,e,bigUp) ...
	     & ~mod(bigS+p+bigUp,2) & ~mod(bigSp+p+bigU,2) ...
	     & ~mod(bigS+e+bigSp,2) & ~mod(bigU+e+bigUp,2);
	% Apply the selection rules
	S=bigS(selx);
	Sp=bigSp(selx);
	U=bigU(selx);
	Up=bigUp(selx);
	
	%    disp(sprintf('The sparsity of this index set is %5.2f%s',...
	%		 100*length(S)/(L+1)^4,'%'))
	
	if length(S)>0
	  % You're now claiming that none of the 3j symbols are zero by virtue of
	  % selection rules - although some may still be zero. Load all at the same
	  % time. 
	  bp=repmat(p,size(S));
	  be=repmat(e,size(S));
	  % Calculate all the zero-bottom row Wigner 3j symbols at once
	  % Work from the preloaded column and element vectors
	  ws=zeroj([S ; Sp ; S ; U],[bp ; bp ; be ; be],[Up ; U; Sp ; Up],...
		   dbL,[],C0,S0);
	  w1234=reshape(ws,length(S),4);
	  w1=(2*S+1) .*w1234(:,1);
	  w2=(2*Sp+1).*w1234(:,2);
	  w3=(2*U+1) .*w1234(:,3);
	  w4=(2*Up+1).*w1234(:,4);
	  % Calculate all the Wigner 6j symbols at once 
	  w6=sixj(S,repmat(e,size(S)),Sp,U,repmat(p,size(S)),Up,...
		  dbL,[],C6,S6)';
	  
	  % Now multiply them all together with the boxcar spectrum and sum
	  % them
	  dsum=dsum+(-1)^(p+e)*(2*e+1)*Be(e==ee,:)*sum(w1.*w2.*w3.*w4.*w6);
	  
	  % Excessive verification at the very end
	  % Don't really ever do this
	  if xver==1
	    disp('Excessive verification of the index arrays')
	    % Make sure that threej with bottom zero returns the same thing      
	    difer(w1-(2*S+1).*zeroj(S,p,Up,dbL,[],C0,S0)')
	    difer(w2-(2*Sp+1).*zeroj(Sp,p,U,dbL,[],C0,S0)')
	    difer(w3-(2*U+1).*zeroj(S,e,Sp,dbL,[],C0,S0)')
	    difer(w4-(2*Up+1).*zeroj(U,e,Up,dbL,[],C0,S0)')
	    % This should be the same as
	    % (no need to initialize SSPUUP since most of the rules are 0)
	    index=0;
	    clear SSPUUP
	    for s=0:L
	      for sp=0:L
		for u=0:L
		  for up=0:L
		    if triangle(s,p,up) & triangle(sp,p,u) ...
			  & triangle(s,e,sp) & triangle(u,e,up) ...
			  & ~mod(s+p+up,2) & ~mod(sp+p+u,2) ...
			  & ~mod(s+e+sp,2) & ~mod(u+e+up,2)
		      index=index+1;
		      % Now satisfy all the conditions etc
		      SSPUUP(index,:)=[s sp u up];
		    end
		  end
		end
	      end
	    end
	    difer(SSPUUP-[S Sp U Up])
	  end % of excessive verification of the indices
	end
      end
      % This is it for a single p, index this appropriately
      GpB(:,p+1)=[dsum./KB.^2]';
      % Check the reduction for p=0
      if xver==1
	disp('Excessive verification of the l=0 term for every degree e')
	if p==0
	  d=repmat(0,1,lTHdo);
	  % This her is for all of them - for the double cap, zeros are
          % in there anyway
	  for e=ee
	    for in=0:L
	      d=d+Be(e==ee,:)*((2*e+1)*(2*in+1)*(2*[0:L]+1)*...
			       [wigner0j(L,e,in).^2]');
	    end
	  end
	  d=[d./KB.^2]';
	  difer(GpB(:,1)-d)
	end
      end
      GpB(:,p+1)=GpB(:,p+1); % This is for a whole set of TH's with meth 1
     case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % This is the Gauss-Legendre integration method
      % Before going in, compute the Gauss-Legendre weights and
      % coefficients since they are reused for every p - good for speed
      [w,x]=gausslegendrecof([2*L+p]);
      GLwc=[w x];
   
      % Initialize this vector once before every degree p
      Gab=repmat(0,Lpot,Lpot);
      % disp(sprintf('Gab initialized for degree p = %i',p))

      % Which degree are we at?
      for mp=-p:p  % Granted, could have switched the order of these
                   % loops, but not if you want conditions like whichb
	% This is the order of the complex spherical harmonic by which we
        % are multiplying both real taper functions... this dictates the
        % loop on a and b... if mp==0 then only a==b. if mp=1 then only 
	% is nonzero
	for a=1:Lpot
	  % Which order are we at with this taper?
	  m1=EMGB(a);
	  % From this and which mp we are at we can figure out over which
          % b's to sum... don't forget this will be symmetric of course
	  if xver~=1
	    if p==0 % FJS & mp==0 
	      whichb=a;
	    else
	      % Next line tries to be smart
	      whichb=find(abs(EMGB)+abs(m1)==abs(mp) | ...
			  abs(abs(EMGB)-abs(m1))==abs(mp));
	      whichb=whichb(whichb>=a)';
	    end
	  else % The excessively loopy method, tries all of them
	    whichb=a:Lpot; 
	  end
	  for b=whichb
	    % Which order are we at with this taper?
	    m2=EMGB(b);
	    % What are we integrating? And remember the REDEFINITION of
            % the normalization
	    % May shorten this - obviously GLMALPHA is loaded MULTIPLE
            % TIMES - TWICE IN ONE RUN, ACTUALLY. Watch with THdo etc.
	    integrand=inline(sprintf(...
		['rindeks(galpha(%i,%i,%i,acos(x),NaN),%i).*'...
		 'rindeks(galpha(%i,%i,%i,acos(x),NaN),%i).*'...
		 'xlm(%i,%i,acos(x))*(-1)^%i*%s'],...
		90*(sord==2)+(-1)^(sord+1)*THdo,L,sord,a,...
		90*(sord==2)+(-1)^(sord+1)*THdo,L,sord,b,...
		p,-mp,mp,newdef));
	    % Now need the longitudinal integrals with the normalization
	    % Most of them are zero in this degenerate case, which is
            % nice. NOW TAKING OUT THE LONGITUDINAL INTEGRAL IS ONLY GOOD
            % FOR THE AXISYMMETRIC CASE, BUT NEEDN'T EVEN ENTER THIS LOOP
            % ON MP IF WE FIGURE OUT WHICH ONES CONTRIBUTE
	    % The following restrictions are less strong than the
            % selection rules; they rather rely on "double-angle formulas"
            % and "trigonometric addition formulas" (see Mathworld)
	    if abs(m1)+abs(m2)==abs(mp) || abs(abs(m1)-abs(m2))==abs(mp)
	      if m1==0 && m2==0
		if mp~=0 
		  absI=0;
		elseif mp==0
		  absI=2*pi;
		end
	      elseif m1==0 || m2==0
		if mp~=0
		  absI=pi;
		elseif mp==0
		  absI=0;
		end
	      elseif mp==0
		if m1*m2>0
		  absI=pi;
		else
		  absI=0;
		end
	      else
		absI=pi/2;
	      end
	    else
	      absI=0;
	    end

	    % Put in the sqrt(2) if needed
	    absI=absI*sqrt(2-[m1==0])*sqrt(2-[m2==0]);
	    % Now add in the colatitudinal integral if the above wasn't zero
	    if absI~=0
	      % This is the colatitudinal integral, in absolute value
	      % Watch the polynomial accuracy of the integration
	      % This makes a HUGE difference if you get it wrong
	      if xver==1
		intgab=absI*abs(gausslegendre([-1 1],integrand,[2*L+p]));
		if p==0 && mp==0 && a==b % Better check this
		  difer(intgab-1/sqrt(4*pi)*eval(newdef))
		end
	      else
		% If mp==0 and a==b then we know the answer, right?
		if p==0 && mp==0 && a==b % Save yourself the trouble
		  intgab=1/sqrt(4*pi)*eval(newdef);
		else
		  % This is when you really do it
		  intgab=absI*abs(gausslegendre([-1 1],integrand,GLwc));
		end
	      end
	    else
	      intgab=0;
	    end
	    
	    if xver==1
		disp(sprintf(...
		    ['a= %2.2i b = %2.2i p = %2.2i m1 = %+i m2 = %+i mp = %+i ',...
		     'intgab = %8.3f'],...
		    a,b,p,m1,m2,mp,(2-[mp==0])*intgab))	      
		if difer(intgab)
		  % Could this still be zero even if absI wasn't, yes,
		  % clearly!! Sometimes the gausslegendre integral is
		  % simply zero. But we should have caught most of them
		  % Remember the tapers (column dimensions) are always
		  % blocksorted even if the degree/order representation
		  % isn't... so this comes out as 0:L, 1:L, -1:L etc for m1
		  % etc 
		  numnz=numnz+1;
		else
		  numz=numz+1;
		  % Which ones after all are still zero?
		  % disp(sprintf(['a= %2.2i b = %2.2i p = %2.2i'],...
		  % ['m1 = %+i m2 = %+i mp = %+i intgab = %8.3f'],...
		  % a,b,p,m1,m2,mp,(2-[mp==0])*intgab))
		end
	    end
	    % For this particular degree p and order mp, here is, we've
            % already taken the absolute value up there - put this in the
            % right place now
	    Gab(a,b)=Gab(a,b)+intgab.^2;
	    % And symmetrize to avoid doing the extra work
	    Gab(b,a)=Gab(a,b);
	  end % Loop over taper index beta
	end % Loop over taper index alpha
	% First I thought only one of them should
        % be nonzero, i.e. that for which mp=-m1-m2... NO that is not
        % true since we are multiplying real functions with a complex
        % harmonic... the Wigner 3j relation doesn't hold in this
        % case. But on the other hand, it is abs(m) which counts
        % here... since the real to complex transformation constructs
        % degree m out of +/- m from the real functions!
      end % Loop over order mp
      % And now we have one element of gammap, but normalize properly!
      Gab=Gab/(2*p+1);

      % So now you have one for GpB - and this should be your answer
      GpB(p+1)=V*Gab*V'/KB^2; % This is for one TH at a time with meth 2
     case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % This is the method using Wigner 3j symbols instead of integration
      Gab=repmat(0,Lpot,Lpot);
      % disp(sprintf('Gab initialized for degree p = %i',p))

      % Define the index arrays for the summation
      % Those are all the possible combinations of l1m1 and l2m2
      % Can use repmat and gamini or meshgrid if you want
      EL1=indeks(repmat(EL(:)',Lpot,1),':');
      EM1=indeks(repmat(EM(:)',Lpot,1),':');
      EL2=repmat(EL(:),Lpot,1);
      EM2=repmat(EM(:),Lpot,1);

      % Now find a way to remember which taper degrees we're using
      [IND1,IND2]=deal(1:length(EL)); % Just a simple index array
      IND1=indeks(repmat(IND1(:)',Lpot,1),':');
      IND2=repmat(IND2(:),Lpot,1);

      % Remember we're working at a single degree p
      % The only non-zero order mp that we require here is M1+M2
      % Why? Check the Gaunt coefficients, it turns into -EMP and
      % subsequently satisfies the order selection rules
      % NOTE, and SEE GAMMAB: THIS IS STRICTLY SPEAKING WRONG, BUT THE
      % END RESULT, WHEN SUMMED OVER THE TAPERS, IS UNCHANGED. HOWEVER,
      % IN GAMMAB WE GET IT RIGHT
      EMP=EM1+EM2; % Really, need the plus here, people

      % The ones that make sense, that is, for this p
      % But of course this needs to be a physical possibility
      qm=abs(EMP)<=p;

      % Apply the selection rules for the zero-bottom coefficients
      no=~mod(p+EL1+EL2,2);
      % And those satisfying the triangle rule for all the others
      tr=triangle(p,EL1,EL2);
      
      % Which of these pass all the checks?
      % These are, for this p, the only sums involved, regardless of the
      % value of alpha, beta, or which tapers applied, etc
      pchk=no & tr & qm;
      
      % Apply all the rules for this case - should move this up
      L1=EL1(pchk);  M1=EM1(pchk);
      L2=EL2(pchk);  M2=EM2(pchk);
      MP=EMP(pchk);  
      ND1=IND1(pchk); ND2=IND2(pchk);
      
      % Without further modification this is always the same
      % So really should be taking this out of the loop

      % These are all the Gaunt coefficients that we will ever need for
      % this degree p - watch out, need to do repmat or result differs!
      GC=gaunt(repmat(p,size(L1)),L1,L2,...
	       MP,M1,M2,[],...
	       dbL,C0,S0,C3,S3)/sqrt(2*p+1);
      % Check that we are indeed picking out the right thing here
      if xver==1
	difer(L1(:)-ELG(ND1(:)))
	difer(M1(:)-EMG(ND1(:)))
	difer(L2(:)-ELG(ND2(:)))
	difer(M2(:)-EMG(ND2(:)))
      end

      % Now we need to multiply in the taper coefficients
      % Once for every taper alpha
      for a=1:Lpot
	% The degrees that matter for this taper are
	% abs(EMGB(a)):L
	% The orders that matter for this taper are
	% [-abs(EMGB(a)) abs(EMGB(a))]
	% Remember if we blocksort, it's for abs(m) now 
	% Once for every taper beta
	whichb=a:Lpot;
	% Should be able to seriously cut this loop here
	for b=whichb
	  % Now all the same with b so we need from the big array
	  % which has all the legal combinations, i.e.
	  % sels=abs(EM1)==abs(EMGB(a)) & abs(EM2)==abs(EMGB(b));
	  % [EL1(sels) EM1(sels) EL2(sels) EM2(sels)]
	  % And this of course is in order
	  % GLMAB(~~GLMAB(:,a),a)
	  % GLMAB(~~GLMAB(:,b),b)
	  % And one repmat and one gamini?

	  % Now we're at an alpha and a beta, and we can do all the
          % coefficients that we have, or: we can also be smarter and a
          % priori exclude the coefficients we know are zero
	  
	  if p==0
	    % If p==0 then all the entries of this guy should be 1/sqrt(4*pi)
	    difer(abs(GC)-1/2/sqrt(pi),[],[],NaN)
	    % Once this is established, can just trust it and not bother
            % to enter this loop
	  end
	  
	  % Just stick in here what you've put in GAMMAB... and test thoroughly

	  % Now perform the actual calculation
	  Gab(a,b)=abs(sum(GLMAB(ND1(:),a).*GC(:).*GLMAB(ND2(:),b))).^2;
	  % Rounding errors no doubt come from sqrt(2) in UMMP
	  % And symmetrize 
	  Gab(b,a)=Gab(a,b);
	  end % Of loop over b
	end % Of loop over a

    % So now you have one for GpB
    GpB(p+1)=V*Gab*V'/KB^2; % This is for a single TH at a time with meth 3
    end % End of switch statement
  end % Loop over degree p i.e. p 

  % Get rid of the NaN's if you calculated only the evens
  GpB=GpB(:,1:evens+1:end);
  % And returns all the p's used, not just the last one
  p=0:evens+1:2*L;

  % End of computation, save them now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % You've done all the THdo at once, but save them separately
  for ind=1:lTHdo
    Gp=GpB(ind,:); 
    K=KB(ind);
    switch isk 
     case 0
      fnpl=sprintf(['%s/%i/GAMMAP-%i-%i-%i-%i-%i.mat'],...
		   fullfile(getenv('IFILES'),'GAMMAP','LTH'),...
		   L,THdo(ind),L,sord,evens,meth);
     case 1
      fnpl=sprintf(['%s/GAMMAP-%i-%i-%i-%i-%i.mat'],...
		   fullfile(getenv('IFILES'),'GAMMAP','LK'),...
		   K,L,sord,evens,meth);
     case 3
      fnpl=sprintf(['%s/GAMMAP-%s-%i-%i-%i.mat'],...
		   fullfile(getenv('IFILES'),'GAMMAP','REGIONS'),...
		   THdo,L,evens,meth);
    end
    % And save this database
    % eval(sprintf('save %s Gp p K',fnpl))
    % Don't save the TH or you're in trouble, later... can always
    % calculate of course
    save(fnpl,'Gp','p','K')
  end

  if xver==1
    disp(sprintf('Number of loops entered that were zero is %i',numz))
    disp(sprintf('Number of loops entered that were nonzero %i', numnz))
  end
  % But you also need to return these guys meshed in with the ones you
  % already knew
  % WATCH OUT - THE RESULT OF INTERSECT WILL BE SORTED
  switch isk
   case 0
    b=ismember(TH,THdo);
   case 1
    b=ismember(TH,KB);
   case 3
    b=1;
  end
  GpBndo(b,:)=GpB;
  KBndo(b,:)=KB;
  Gp=GpBndo; 
  K=KBndo;
else
  % Return whatever you got if you don't need to do any more
  Gp=GpB;
  K=KB;
end

% Prepare output
varns={Gp,p,K};
varargout=varns(1:nargout);

elseif strcmp(L,'demo1')
  % Compare SINGLE CAP results for different methods - SLOW!!!
  Gp1=gammap(3,65,1,0,1); 
  Gp2=gammap(3,65,1,0,2);
  Gp3=gammap(3,65,1,0,3);
  difer(Gp1-Gp2)
  difer(Gp1-Gp3)
  difer(Gp2-Gp3)
elseif strcmp(L,'demo2')
  % Compare double cap results for different methods - SLOW!!!
  Gp1=gammap(3,65,2,0,1); 
  Gp2=gammap(3,65,2,0,2);
  Gp3=gammap(3,65,2,0,3); 
  difer(Gp1-Gp2)
  difer(Gp1-Gp3)
  difer(Gp2-Gp3)
  % Check that single and double caps are consistent
elseif strcmp(L,'demo3')
  difer(gammap(10,180,1)-gammap(10,0,2))
elseif strcmp(L,'demo4')
  L=29;
  Gp1=gammap(4,L,1,0,1);
  Gp2=gammap(4,L,1,0,2);
  Gp3=gammap(4,L,1,0,3);
  difer(Gp1-Gp2)
  difer(Gp1-Gp3)
  difer(Gp2-Gp3)
  Gp1=gammap(4,L,2,0,1);
  Gp2=gammap(4,L,2,0,2);
  Gp3=gammap(4,L,2,0,3);
  difer(Gp1-Gp2)
  difer(Gp1-Gp3)
  difer(Gp2-Gp3)
end

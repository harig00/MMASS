function varargout=gammab(L,TH,sord,evens,meth,isk)
% [Gab,p,K,V]=GAMMAB(L,TH,sord,evens,meth,isk)
%
% Like GAMMAP except it retains the individual (alpha,beta) pairs.
% Dahlen and Simons 2008 eq. (164)
%
% INPUT:
%
% L         Taper bandwidth that is required in this case, always a scalar
% TH        Colatitudinal radius of the cap, in degrees <=180, may be vector
%           Colatitudinal halfwidth of the cut, degrees <=90, may be vector
%           ... or: the Shannon number... see option 'isk'
% sord      1 Single cap of diameter 2TH [default]
%           2 Double cap left by subtracting belt of width 2TH (watch out!)
% evens     1 only the even p's [default, good for variance but not covariance]
%           0 all p's 
% meth      1 Not available
%           2 using GL integration [exploiting symmetry of assumed caps]
%           3 method using Wigner 3j symbols [adaptable to irregular regions]
% isk       0 nothing unusual here [default]
%           1 the TH is actually a Shannon number K, and TH is calculated
% 
% OUTPUT:
%
% Gab       Gamma(alpha,beta,p) as defined by Dahlen & Simons (2008)
%           Dimensions are: length(TH) x (L+1)^2 x (L+1)^2 x length(p)
% p         The degrees at which this is evaluated, containing evens or not
% K         The Shannon number(s) belonging to the individual domains
% V         The eigenvalues belonging to the individual domains
%
% EXAMPLE & VERIFICATION
% 
% GAMMAB('demo1',TH,sord,evens) % GAMMAB vs GAMMAP, meth 2 vs meth3
%
% Last modified by fjsimons-at-alum.mit.edu, 04/25/2007

tic
defval('L',10)
defval('TH',30)
defval('sord',1)
defval('evens',0)
defval('meth',3)
defval('isk',0)

if ~isstr(L)

if meth==1  
 error('I just left the option method 1 in but only 2 and 3 are allowed')
end

% The dimension of the basis
Lpot=(L+1)^2;

% Initialize four-dimensional array
GabB=repmat(NaN,[length(TH) Lpot Lpot 2*L+1]);
GabB=GabB(:,:,:,1:evens+1:end);
KB=repmat(NaN,length(TH),1);
VB=repmat(NaN,length(TH),Lpot);

% For multiple integer THs
ndto=0;
for ind=1:length(TH)
  switch isk
    case 0
     fnpl=sprintf(['%s/%i/GAMMAB-%i-%i-%i-%i-%i.mat'],...
		  fullfile(getenv('IFILES'),'GAMMAB','LTH'),...
		  L,TH(ind),L,sord,evens,meth);
     % Perhaps you wanted only the evens and you had them all?
     fnplb=sprintf(['%s/%i/GAMMAB-%i-%i-%i-%i-%i.mat'],...
		  fullfile(getenv('IFILES'),'GAMMAB','LTH'),...
		  L,TH(ind),L,sord,evens-1,meth);
     case 1
     fnpl=sprintf(['%s/GAMMAB-%i-%i-%i-%i-%i.mat'],...
		  fullfile(getenv('IFILES'),'GAMMAB','LK'),...
		  TH(ind),L,sord,evens,meth);
   otherwise
    error('Specify valid option for isk')
  end
  if exist(fnpl,'file')==2
    % Got these
    load(fnpl)
    disp(sprintf('Loading %s',fnpl))
    GabB(ind,:,:,:)=Gab;
    KB(ind,:)=K;
    VB(ind,:)=V;
  elseif exist(fnplb,'file')==2
    load(fnplb)
    GabB(ind,:,:,:)=Gab(:,:,1:evens+1:end);
    p=p(1:evens+1:end);
    KB(ind,:)=K;
    VB(ind,:)=V;
  else
    disp(sprintf('%s does not exist',fnpl))
    % Still need these
    ndto=ndto+1;
    THdo(ndto)=TH(ind);
  end
end

% Still need to do a bunch here; do them all at once
if ndto>=1
  if ndto>1 && meth~=1
    error(sprintf('Method %i only for (additional) scalar TH',meth))
  end
  % Don't do these - the missing ones have NaN's in them
  GabBndo=GabB;
  KBndo=KB;
  VBndo=VB;

  switch isk
   case 0
    % Calculate the area from the input angle TH
    A=4*pi*spharea(THdo,sord);
    if sord==3
      error('This choice is not allowed')
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
      error('This choice is not allowed')
    end
    % Watch out for formally correct but meaningless results
    KB=KB(~imag(THdo) & THdo>=0);
    THdo=THdo(~imag(THdo) & THdo>=0);
    % If it doesn't make sense, the code returns NaN and only the sensible ones
    % will be saved. What if both empty? Do nothing
    if isempty(THdo) && isempty(KB)
      % Return whatever you got if you don't need to do any more
      Gab=GabB;
      K=KB;
      return
    end
  end
  
  % Re-initialize after possible deletion of certain imaginary THdo as
  % before - but remember, length(THdo) is only allowed to be 1 here
  GabB=repmat(NaN,[Lpot Lpot 2*L+1]);

  % Get the power spectrum of the boxcar cap
  [Be,ee]=bpboxcap(THdo,2*L,[],sord-1,sord);
  % If it's a DOUBLE cap, will need/get only even degrees e
  
  % What are the combinations of 4 indices?
  bigS=        gamini([0:L],(L+1)^3)';
  bigSp=repmat(gamini([0:L],(L+1)^2)',(L+1)^1,1);
  bigU= repmat(gamini([0:L], L+1)',   (L+1)^2,1);
  bigUp=repmat(       [0:L]',         (L+1)^3,1);
  
  if meth~=1
    % Get the coefficients for the tapers - you always need the
    % eigenvalues
    % Watch the conventions of GLMALPHA and GRUNBAUM2 have different conventions
    blox=0;
    [GLMAB,VB,ELG,EMG]=glmalpha(90*(sord==2)+(-1)^(sord+1)*THdo,L,sord,blox);
    newdef='4*pi';
    disp(sprintf('GAMMAP Redefined normalization for tapers to %s',newdef))
    GLMAB=GLMAB*sqrt(eval(newdef));
    % Get some other things too
    [EM,EL,mz,blkm]=addmout(L);
    % Block sorted array for future reference - could have pulled it out
    % of GLMALPHA since its columns are always sorted this way for
    % axisymmetric regions
    EMGB=EM(blkm); % Only used for meth 2 anyway 
  end

  if meth~=2
    % Must load EXACTLY the one with the specified bandwidth, no default
    % substitutions are to be accepted, if not, won't know which vector
    % you've got!!! So this is how to do this
    dbL=2*L;
    % Load the ZEROJ database just once, for the maximum at all required
    [jk,C0,S0]=zeroj(0,0,0,dbL);
    % Load the THREEJ database just once, as well
    [jk,C3,S3]=threej(0,0,0,0,0,0,dbL);
  end

  if meth==3
    % Don't forget to make them complex for the Wigner 3j formalism to work
    U=ummp(L);
    % Note that you now give up the strict identification with +/- m
    % It's indeed destroying the block-sorting, rather now going by
    % larger blocks of +/- m
    GLMAB=U*GLMAB;
  end

  disp(sprintf('GAMMAB Using METHOD %i',meth))
  
  % There's something funny - the first time round, if the database of 0j
  % or 3j isn't there, it makes it, but then it rereads it at every step
  % in the loop rather than using the loaded database... for some reason

  % Initialize this matrix here so you can save it for later
  Gab=repmat(0,[Lpot Lpot 2*L+1]);
  
  % For every desired degree p
  % Make p=0 a special case but don't adapt the inside of the loop which
  % would confirm it. This is special to GAMMAB
  p=0;
  disp(sprintf('Working on %i of %i',p,2*L))
  Gab(:,:,p+1)=4*pi*eye(Lpot,Lpot);
  % Go on for higher p
  for p=evens+1:evens+1:2*L
    disp(sprintf('Working on %i of %i',p,2*L))
    switch meth
     case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      error('Method 1 is not possible here')
     case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % This is the Gauss-Legendre integration method
      % Before going in, compute the Gauss-Legendre weights and
      % coefficients since they are reused for every p 
      [w,x]=gausslegendrecof([2*L+p]);
      GLwc=[w x];
      
      for q=-p:p 
      % Which degree are we at?
      for a=1:Lpot
	% Which order are we at with this taper? This assumes symmetry
	m1=EMGB(a);
	% Which b's do we need?  Orthogonality
	if p==0 % FJS & q==0
	  whichb=a;
	else
	  % Next line tries to be smart but is caught later anyway
	  whichb=find(abs(EMGB)+abs(m1)==abs(q) | ...
		      abs(abs(EMGB)-abs(m1))==abs(q));
	  whichb=whichb(whichb>=a)';
	end
	for b=whichb
	    % Which order are we at with this taper?
	    m2=EMGB(b);
	    % Get rid of this - this repeated loading is ridiculously slow
	    integrand=inline(sprintf(...
		['rindeks(galpha(%i,%i,%i,acos(x),NaN),%i).*'...
		 'rindeks(galpha(%i,%i,%i,acos(x),NaN),%i).*'...
		 'xlm(%i,%i,acos(x))*(-1)^%i*%s'],...
		90*(sord==2)+(-1)^(sord+1)*THdo,L,sord,a,...
		90*(sord==2)+(-1)^(sord+1)*THdo,L,sord,b,...
		p,-q,q,newdef));
	    % Now need the longitudinal integrals with the normalization
	    if abs(m1)+abs(m2)==abs(q) || abs(abs(m1)-abs(m2))==abs(q)
	      if m1==0 && m2==0
		if q~=0 
		  absI=0;
		elseif q==0
		  absI=2*pi;
		end
	      elseif m1==0 || m2==0
		if q~=0
		  absI=pi;
		elseif q==0
		  absI=0;
		end
	      elseif q==0
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
	    % [p q m1 m2 absI]
	    
	    % Now add in the colatitudinal integral if the above wasn't zero
	    if absI~=0
	      % This is the colatitudinal integral, in absolute value
	      % If q==0 and a==b then we know the answer, right?
	      if p==0 && q==0 && a==b % Save yourself the trouble
		intgab=1/sqrt(4*pi)*eval(newdef);
	      else
		intgab=absI*abs(gausslegendre([-1 1],integrand,GLwc));
	      end
	    else
	      intgab=0;
	    end
	    
	    % For this particular degree p and order q, here it is
	    Gab(a,b,p+1)=Gab(a,b,p+1)+intgab.^2;
	    % And symmetrize to avoid doing the extra work - use tril?
	    Gab(b,a,p+1)=Gab(a,b,p+1);
	  end % Loop over taper index beta
	end % Loop over taper index alpha
      end % Loop over order q
      % And now we have one element of gammap, but normalize properly!
      Gab(:,:,p+1)=Gab(:,:,p+1)/(2*p+1);

     case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Define the index arrays for the summation of both sets of tapers
      EL1=indeks(repmat(EL(:)',Lpot,1),':');
      EM1=indeks(repmat(EM(:)',Lpot,1),':');
      EL2=repmat(EL(:),Lpot,1);
      EM2=repmat(EM(:),Lpot,1);

      % Now find a way to remember which taper indices we need
      [IND1,IND2]=deal(1:length(EL)); % Just a simple index array
      IND1=indeks(repmat(IND1(:)',Lpot,1),':');
      IND2=repmat(IND2(:),Lpot,1);

      % Remember we're working at a single degree p
      EQ=EM1+EM2; % Really, need the plus here, people due to GAUNT setup
      % Thee complex taper coefficients have lost their strict
      % identification with the orders, rather they are now identified
      % with the ABSOLUTE order. 
      % The ones that make sense, that is, for this p
      qm=abs(EQ)<=p;

      % Apply the selection rules for the zero-bottom coefficients
      no=~mod(p+EL1+EL2,2);
      % And those satisfying the triangle rule for all the others
      tr=triangle(p,EL1,EL2);
      
      % Which of these pass all the checks?
      pchk=no & tr & qm;
      
      % Apply all the rules for this case
      L1=EL1(pchk);  M1=EM1(pchk);
      L2=EL2(pchk);  M2=EM2(pchk);
      Q=EQ(pchk);  
	  
      % Maybe this next line is wrong, or else the sign of the GC?
      ND1=IND1(pchk); ND2=IND2(pchk);

      % Check that we are indeed picking out the right thing here
      % difer(L1(:)-ELG(ND1(:)))
      % difer(M1(:)-EMG(ND1(:)))
      % difer(L2(:)-ELG(ND2(:)))
      % difer(M2(:)-EMG(ND2(:)))

      % These are all the Gaunt coefficients that we will ever need
      GC=gaunt(repmat(p,size(L1)),L1,L2,...
	       Q,M1,M2,[],...
	       dbL,C0,S0,C3,S3)/sqrt(2*p+1);
      
      % Now we need to multiply in the taper coefficients
      for a=1:Lpot
	% Once for every taper beta except when p==0
	if p~=0
	  whichb=a:Lpot;
	else
	  whichb=a;
	end
	for b=whichb
	  % CHECK GAUNT/THREEJ POLICY ON REPEATED INDICES
	  % Hold on - looks like this is sensitive to THE ORDER! and 
	  % Now perform the actual calculation
	  % Look at the common coefficients: they occur for equal abs(m)
	  % [M1 M2 ~~abs(GLMAB(ND1(:),a).*GLMAB(ND2(:),b))]
	  % [repmat(p,size(L1)),L1,L2,Q,M1,M2 ~~abs(GLMAB(ND1(:),a).*GLMAB(ND2(:),b))]
	  % In the sum, some are over m's and some are over q's
	  % [Q GLMAB(ND1(:),a).*GC(:).*GLMAB(ND2(:),b)]
	  % Need to group them by Q, sum within them, abs, ^2, sum again
	  % Get an intermediate result
	  intres=GLMAB(ND1(:),a).*GC(:).*GLMAB(ND2(:),b);
	  % Figure out which ones are a priori interesting - i.e. nonzero
	  ofintrest=~~abs(intres);
	  % For a general area, one would not expect many of these to be
          % nonzero-the sum would a priori involve all orders. 
	  % But I don't want to go any further here in restricting this
          % loop to ANTICIPATE which ones are going to be zero - method 3
          % will be the one used for arbitrarily shaped areas.
	  % And sum the magnitude of this per q, e.g. for axisymmetric
          % regions when only a single when abs(m)=1 is
          % involved, we get contributions from q=-2, 2, and 0
	  if ~~sum(ofintrest)
	    Qinv=Q(ofintrest);
	    intres=intres(ofintrest);
	    % The next line is probably overkill of sophistication
	    % Find the ranges belonging to same q
	    okill=sortrows([Qinv intres]);
	    [ol,bol,snol]=degamini(okill(:,1));
	    % Use sparse array to do the blocky stuff... really should turn
	    % this into a function to compare with ROW2STATS and BLOCKMEAN
	    % But this should be the right result
	    Gab(a,b,p+1)=sum(abs(sparse(gamini(1:length(ol),bol),...
					1:snol(end),1)*okill(:,2)).^2);
	    % The fact that I was doing this wrong internally inside
            % GAMMAP had no detrimental effect on the end result.
	    % And symmetrize - perhaps should be doing this at the very end?
	    Gab(b,a,p+1)=Gab(a,b,p+1);
	  end
	  end % Of loop over b
	end % Of loop over a
    end % End of switch statement
  end % Loop over degree p i.e. p 

  % And returns all the p's used, not just the last one
  p=0:evens+1:2*L;

  % End of computation, save them now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Get rid of the NaN's if you calculated only the evens
  GabB=Gab(:,:,1:evens+1:end); % New compared to GAMMAP

  % In this code length(THdo) can only be one so there is no more loop
  Gab=GabB(:,:,:); % Straight from the calculation
  K=KB;
  V=VB;
  switch isk 
   case 0
    fnpl=sprintf(['%s/%i/GAMMAB-%i-%i-%i-%i-%i.mat'],...
		   fullfile(getenv('IFILES'),'GAMMAB','LTH'),...
		 L,THdo,L,sord,evens,meth);
   case 1
    fnpl=sprintf(['%s/GAMMAB-%i-%i-%i-%i-%i.mat'],...
		 fullfile(getenv('IFILES'),'GAMMAB','LK'),...
		 K,L,sord,evens,meth);
  end
  
  % Before you save this, set the tiny ones to really zero
  Gab(Gab<eps)=0;
  GabB(GabB<eps)=0;
  
  keyboard
  
  % And save this database
  save(fnpl,'Gab','p','K','V')
  
  switch isk
   case 0
    b=ismember(TH,THdo);
   case 1
    b=ismember(TH,KB);
  end
  GabBndo(b,:,:,:)=GabB;
  KBndo(b,:)=KB;
  VBndo(b,:)=VB;
  Gab=GabBndo; 
  K=KBndo;
  V=VBndo;
else
  % Return whatever you got if you don't need to do any more
  Gab=GabB;
  K=KB;
  V=VB;
end

varns={Gab,p,K,V};
varargout=varns(1:nargout);

% % ALWAYS check that the eigenvalue-weighted sum returns GAMMAP!
if L<=20
  [Gp,p,K]=gammap(L,TH,sord,evens,1);
else
  [Gp,p,K]=gammap(L,TH,sord,evens,3);
end
for int=1:length(TH)
  for ind=1:length(p)
    difer(V(int,:)*squeeze(Gab(int,:,:,ind))*V(int,:)'...
	  /K(int,:)^2-Gp(int,ind),[],[],NaN)
  end
end
disp('GAMMAB Checked against GAMMAP')

elseif strcmp(L,'demo1')  
  L=4;
  [Gab3,p3,K3,V3]=gammab(L,TH,sord,evens,3);
  Gab3=squeeze(Gab3);
  [Gab2,p2,K2,V2]=gammab(L,TH,sord,evens,2);
  Gab2=squeeze(Gab2); 
  % Compare with method 1 for GAMMAP which is independent and fast
  [Gp,p,K]=gammap(L,TH,sord,evens,1);
  for ind=1:length(p)
    difer(V2*Gab2(:,:,ind)*V2'/K2^2-Gp(ind),[],[],NaN)
    difer(V3*Gab3(:,:,ind)*V3'/K3^2-Gp(ind),[],[],NaN)
  end
  % Also check the results among themselves
  difer(Gab2-Gab3)
end

% TO DO:
% COMPARE THE LAST TWO OR THREE THINGS I HAVE OF METHOD 3 WITH METHOD 2
% CHANGE THE LOADING OF THE GALPHA'S IN METHOD 2... FOR REAL!
% AFTHER THAT, WE HAVE ONE MEMORY AND ONE SPEED-EFFICIENT WAY FOR THE
% COMPUTATION - THE MEMORY-EFFICIENT WAY EXPLOITS THE SYMMETRY AND CAN'T
% BE ADAPTED TO IRREGULAR REGIONS. MAYBE FURTHER EXPLOITATION IS ON THE
% WAY. CHECK HOW MANY REPEATED ENTRIES EVERY SUBMATRIX HAS!

toc

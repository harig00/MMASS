function [Cab,morder,N]=whitevariance(l,L,TH,meth,memeff)
% [Cab,morder]=WHITEVARIANCE(l,L,TH,meth,memeff)
%
% Computes the multitaper variance for white source processes and
% axisymmetric Slepian functions, at a SINGLE DEGREE. Single polar cap.
% Dahlen & Simons (2008) eq. (162-163)
% 
% INPUT:
%
% l        The degree at which you want this evaluated
% L        The bandwidth of the tapers used
% TH       The colatitudinal radius of the polar cap
% meth     0 Using Wigner-3j coefficients computed on the fly
%          1 Using Wigner-3j coefficients previously computed [default]
%          2 Using Gaunt coefficients and Gauss-Legendre integration
%          0->1 can be adapted for arbitrarily shaped areas
%          2   is only for axisymmetric areas (degenerate orders) but can
%          be made even faster by not entering a loop on mp 
% memeff   Only for meth=1: 0 Elegant formulation that is memory-intensive
%                           1 Memory-efficient version [default]
%
% OUTPUT:
%
% Cab      The spherical multitaper (co)variance matrix
% morder   The order of the tapers in the order they are presented
% N        The Shannon number
%
% EXAMPLE:
%
% difer(whitevariance(0,1,30,1)-whitevariance(0,1,30,0))
% difer(whitevariance(0,1,30,1)-whitevariance(0,1,30,2))
% difer(whitevariance(2,2,30,1)-whitevariance(2,2,30,0))
% difer(whitevariance(2,2,30,1)-whitevariance(2,2,30,2))
% difer(whitevariance(4,3,30,1)-whitevariance(4,3,30,0))
% difer(whitevariance(4,3,30,1)-whitevariance(4,3,30,2))
%
% difer(whitevariance(15,10,20,1)-multivar(10,20,1,15))
% difer(whitevariance(40,3,30,1)-multivar(3,30,1,40));
% difer(whitevariance(25,24,30,1)-multivar(24,30,1,25))
%
% Last modified by fjsimons-at-alum.mit.edu, 04/29/2007
%
% This works, although inefficiently - compare with GAMMAB and MULTIVAR

% MULTIVAR had a special problem to solve since it disentangled the
% matrix GAMMAB from the expression and needed to account for this using
% the fixed trick with the sparse matrix... but both approaches are
% equivalent.

% May load databases and pass them on? Loading takes time.

defval('l',2)
defval('L',3)
defval('TH',30)
defval('meth',1)
defval('memeff',1)

% Excessive testing? A good idea, usually;
defval('xver',0)

fnpl=sprintf('%s/whitevariance-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'VARIANCE'),l,L,TH,meth);

if exist(fnpl,'file')==2 
  disp(sprintf('load %s',fnpl))
  eval(sprintf('load %s',fnpl))
else
  t=cputime;

  % Axisymmetric polar cap only, for the time being
  % Get the localization functions; sort as in ADDMOUT
  blox=0;
  sord=1;
  [GLMAB,V,ELG,EMG,N]=glmalpha(TH,L,sord,blox);
  % SHOULD REALLY RENORMALIZE THEM HERE AS IN GAMMAB
  % Make these coefficients complex now
  U=ummp(L);
  % Do it!
  GLMAB=U*GLMAB;
  
  Lpot=(L+1)^2;
  [EM,EL,mz,blkm]=addmout(L);
  % Block sorted array for future reference 
  ELGB=ELG(blkm);
  EMGB=EMG(blkm);
  % Check they are sorted the same as those of GLMALPHA - necessary?
  difer(ELG-EL); difer(EMG-EM); 

  switch meth
   case 1
    % Get orders and degrees not in block form
    % Perhaps put in here min(2*l,2*L)
    [EMP,ELP]=addmout(min(2*l,2*L),1);
    % Specify all degrees and orders of Gaunt coefficients you need...
    % ... for the summation over lpmp (only evens)... already so
    % EMP=EMP(~mod(ELP,2));
    % ELP=ELP(~mod(ELP,2));
    NLP=length(ELP);

    % Define the summation arrays without the obvious zeros
    switch memeff
     case 1
      % Some sort of initialization is required here
      % The factor by which the above memory-inefficient procedure would
      % have wasted the memory of the arrays
      fax=0.0035;  % Taken from a run that broke afterwards
      % May reduce this factor; but if you reduce it too much,
      % the process will run slower; keep track of this by writing it out
      [LP,MP,L1,M1,L2,M2,ind1,ind2]=...
	  deal(repmat(NaN,max(2,round(Lpot^2*NLP*fax)),1));
      % Must put in max 2 if not get wrong dimension!
      runs=0;
      % Looks like we'll be able to take this down a notch - 
      % the lp is also limited not only by 2*l but also by the maximum
      % value that l1 and l2 can obtain, which is 2*L;
      % New thing here... is the min 2*L inclusion
      % But actually, later on, you'll need them... or do you?
      for lp=0:2:min(2*l,2*L)
	for mp=-lp:lp
	  % ACTUALLY, FOR THE SYMMETRIC CASE, SHOULD NOT ENTER INTO A M1
          % M2 LOOP AT ALL. FOR AXISYMMETRIC CASE, SUM IS NOT DOUBLE BUT
          % SINGLE; FOR AXISYMMETRIC CASE, THERE ISN'T EVEN A SUM OVER
          % THE ORDERS M1 and M2
	  % But leave this as it is more general
	  % In other words, this is just based on the selection rules,
          % but we do know in advance that we're dealing with only a
          % handful of non-zero taper coefficients!
	  for one=1:length(EMG)
	    l1=ELG(one);
	    m1=EMG(one);
	    % noselrul=1:length(EMG);
	    selrul=find(~mod(ELG+lp+l1,2) & ...
		EMG==mp-m1 & ELG<=lp+l1 & ELG>=abs(lp-l1));
	    % So the length of selrul is what repeats one
	    ind1(runs+1:runs+length(selrul))=one;
	    % Same convention as down here with the same 
	    ind2(runs+1:runs+length(selrul))=selrul;
	    % This is important
	    for two=selrul(:)'
	      l2=ELG(two);
	      m2=EMG(two);
	      
	      % Collect the arrays that will finally be used
	      runs=runs+1;

	      LP(runs)=lp;
	      MP(runs)=mp;
	      L1(runs)=l1;
	      M1(runs)=m1;
	      L2(runs)=l2;
	      M2(runs)=m2;
	    end
	  end
	end
      end
      disp('Finished triple loop building index matrices')
      L1=L1(~isnan(L1));
      L2=L2(~isnan(L2));
      LP=LP(~isnan(LP));
      M1=M1(~isnan(M1));
      M2=M2(~isnan(M2));
      MP=MP(~isnan(MP));
      ind1=ind1(~isnan(ind1));
      ind2=ind2(~isnan(ind2));
      fax2=length(L1)/(Lpot^2*NLP);
      disp(sprintf('Percentage retained a priori %6.3f',fax*100))
      disp(sprintf('Percentage retained for summation %6.3f',fax2*100))
      disp('Look to improve this ratio from the start')
      fid=fopen(fullfile(getenv('IFILES'),'VARIANCE','wvfactors'),'a+');
      fprintf(fid,'%i %i %8.6f\n',[l L fax2]);
      fclose(fid);
      
      % This is to get NLPR without having the original long sequence
      NLPR=zeros(1,NLP);
      [jk1,alt,jk2]=degamini(MP);
      NLPR(1:length(alt))=alt;
     case 0 
       % ... for the summation over l1m1 and l2m2
       % you may use MESHGRID if you want
       % The size of these arrays will then be (L+1)^4
       EL1=indeks(repmat(EL(:)',Lpot,1),':');
       EM1=indeks(repmat(EM(:)',Lpot,1),':');
       
       EL2=repmat(EL(:),Lpot,1);
       EM2=repmat(EM(:),Lpot,1);
       
       % All together thus: the size is (L+1)^4*NLP
       % Watch out that you do not exceed the maximum variable size
       % Find something better than this
       L1=repmat(EL1(:),NLP,1);
       M1=repmat(EM1(:),NLP,1);
       L2=repmat(EL2(:),NLP,1);
       M2=repmat(EM2(:),NLP,1);
       
       LP=indeks(repmat(ELP(:)',length(EL1),1),':');
       MP=indeks(repmat(EMP(:)',length(EL1),1),':');
       
       % For l1m1, pick kmp out of a sequence in which each of the Lpot
       % degrees individually is repeated Lpot times & find the original index
       % For l2m2, pick kmp out of a Lpot times repeated regular order sequence
       % of Lpot entries & find which original indices are being picked
       
       % These are the original lengths before the selection rules etc
       % disp(sprintf('%i %i %i',length(L1),length(L2),length(LP)))
       
       % Check out the summation arrays:
       % imagesc([LP MP L1 M1 L2 M2])
       % Take only the ones whose degree sums are even
       no=~mod(LP+L1+L2,2);
       % And those satisfying the triangle rule
       tr=triangle(LP,L1,L2);
       % And those satisfying the zero-sum rule
       % WATCH and correct this
       zs=~[-MP+M1+M2];
       kpm=[no & tr & zs];
       clear no zs tr
       
       LP=LP(kpm); MP=MP(kpm);
       L1=L1(kpm); M1=M1(kpm);
       L2=L2(kpm); M2=M2(kpm);
       
       % This is after the selection rules etc
       % disp(sprintf('%i %i %i',length(L1),length(L2),length(LP)))
       
       % This is how many you pick for l2m2 per lpmp
       lmsel=reshape(kpm,Lpot^2,NLP);
       % This is the number of elements for a single lpmp
       % ie the number of times a single lpmp is repeated
       NLPR=sum(lmsel,1);
       difer(sum(NLPR)-sum(kpm)) % Duh
       flmsel=find(lmsel);
           
       % The index of the chosen L1M1 into the original ELGEMG per LPMP 
       ind1=mod(flmsel,Lpot^2); ind1(~ind1)=Lpot^2;
       ind1=ceil(ind1/Lpot);
       
       % The index of the chosen L2M2 into the original ELGEMG per LPMP
       ind2=mod(flmsel,Lpot); ind2(~ind2)=Lpot;
     otherwise
      error('Specify correct memory efficiency parameter')
    end
    difer(gamini(ELP,NLPR)-LP') % Doh

    if xver==1
      % This may only be the case for when p=0?
      difer(ind1(1:Lpot)-[1:Lpot]')
      difer(L1(:)-ELG(ind1(:)))
      difer(M1(:)-EMG(ind1(:)))
      difer(L2(:)-ELG(ind2(:)))
      difer(M2(:)-EMG(ind2(:)))
    end
    % And get ALL of the nonzero Gaunt coefficients at once
    % In the trials just one more of them was zero without "obvious"
    % reason... It was (2,3,3,0,2,-2) => see Edmonds' book
    GC=gaunt(LP,L1,L2,MP,M1,M2,'table')./sqrt(2*LP(:)'+1);
    % Remember the factors sqrt(2l+1) etc are already in there

    disp('Got all Gaunt coefficients of the first part')
    
    % [LP MP L1 ELG(ind1(:)) M1 EMG(ind1(:)) ...
    % L2 ELG(ind2(:)) M2 EMG(ind2(:)) ind2(:) round(GC(:)*1000)]
    % [LP MP L1 M1 L2 M2 round(GC(:)*1000)]
    % and take this nplmpm at a time
    % Make a sparse matrix for the purposes of the summation 
    % And then need to average this over the orders available per lpmp
    % Like BLOCKMEAN... should perhaps make a RANGEMEAN function... but see
    % also ROW2STATS, obviously
    ro=gamini(1:NLP,NLPR);
    co=1:length(GC);
    % Specify size of the matrix in case there were zeroes
    CT=sparse(ro,co,1,NLP,length(GC));
    
    % This is just about the integral over the product of both functions
    % Make a loop over the elements of Cab; note that b needs to be
    % reordered; every alpha or beta belongs to a very specific order, and
    % by picking out its elements ind2 in some order, you are reshuffling
    % this; you thus need to map b to reflect this... NO: it's blocks that
    % are reshuffled, these blocks are L-abs(m)+1 long
    % Initialization saves time
    intg=repmat(0,[Lpot Lpot NLP]);

    for a=1:Lpot
      for b=a:Lpot
	% What's bth in the sequence belongs elsewhere in the matrix
	% Now try to modify this intelligently in case GLMAB is an empty
        % thing due to the symmetry of the polar-cap kernels
	bprime=addmabout(L,ELGB(b),EMGB(b));
	intg(a,b,:)=CT*[GLMAB(ind1(:),a).*GC(:).*GLMAB(ind2(:),bprime)];
      end
    end
    disp('Finished alpha-beta matrix')
    if xver==1
      for index=1:size(intg,3)
	if sum(sum(intg(a,b,index)))==0
	  % Maybe some of this has to do with the degeneracy?
	  disp(sprintf('You''re doing too much work on %i or l''=%i',...
		       index,ELP(index)))
	end
      end
      % If lp=mp=0 the intg matrix should be diagonal! Check this out.
      % Note that the sign of the eigenfunction was chosen arbitrarily and
      % that the inner product may be plus or minus depending on the
      % situation. Note that ylm(0,0,0)=1/sqrt(4*pi) already
      difer(diag(abs(intg(:,:,1))*sqrt(4*pi))-1) % Should be one on the diagonal
      difer(intg(:,:,1)-diag(diag(intg(:,:,1)))) % Should be zero off the diagonal
    end

    % See if you ever need all of them or if you can get by with half of them.
    % If ELP is too big, do it over.... does it mean that the last 3d
    % panels of intg are always zero even if GC isn't? I think so.
    
    % Now the third dimension contains the lpmp... get the Gaunt coefficients
    GC=gaunt(l,l,ELP,0,0,0,'table').*sqrt(2*ELP(:)'+1)/(2*l+1);
    disp('Got all Gaunt coefficients of the second part')

    Cab=repmat(0,Lpot,Lpot);
    for index=1:NLP
      Cab=Cab+GC(index)*abs(intg(:,:,index)).^2;
    end

    % Normalize, let the signal strength be S=1
    Cab=Cab/sqrt(pi);
   case 2 
    % By Gauss-Legendre integration
    % Single or double cap
    sord=1;
    Cab=repmat(0,Lpot,Lpot);
    ELP=0:2:min(2*l,2*L);
    % Get all the zero-bottom Gaunt coefficients at once
    GC=gaunt(l,l,ELP,0,0,0,'table')./sqrt(2*ELP+1)/(2*l+1);
    
    for a=1:Lpot
      for b=a:Lpot
	%      disp(sprintf('Tapers %i %i',a,b))
	%      fid=fopen(sprintf('Taper%i_%i',a,b),'a+');
	% What is the order of these functions? Remember, they were
	% block sorted in the column dimension, regardless of whether
	% they'd been block sorted in the row dimension.  This is the
        % part that only works for axisymmetric functions
	m1=EMGB(a);
	m2=EMGB(b);
	%      disp(sprintf('m1 %+i m2 %+i',m1,m2))
	% So check that the only non-zero coefficients are indeed within
	% that order
	for index=1:length(ELP)
	  lp=ELP(index);
	  % Can I obviate this loop here by going only by the selection
          % rules or what not? These were derived from the longitudinal integral
	  % disp(' ')
	  % See if we get them all with this improvement
	  mpcand=unique([abs(m1)+abs(m2) abs(abs(m1)-abs(m2))]);
	  mpcand=mpcand(abs(mpcand)<=lp);
	  % Don't do negatives, just do these later (look for mp==0)
	  % mpcand=unique([-mpcand mpcand]);
	  % There is a difference between [] and empty
	  if isempty(mpcand); mpcand=[]; end
	  % mpcand=-lp:lp;
          % Also realize that abs(mp) counts, +/- values are identical
	  for mp=mpcand
	    % It's obviously not smart to reload GALPHA all the time, 
	    % but leave it be as it is fast enough
	    % LATER: SWITCH FORMULATION TO LOAD TWO SIMULTANEOUSLY
	    % USING THE COMPLEX CONJUGATE MAKES NO DIFFERENCE FOR YLM
	    integrand=inline(sprintf(...
		['rindeks(galpha(%i,%i,%i,acos(x),NaN),%i).*'...
		 'rindeks(galpha(%i,%i,%i,acos(x),NaN),%i).*'...
		 'xlm(%i,%i,acos(x))*(-1)^%i'],...
		TH,L,sord,a,TH,L,sord,b,lp,-mp,mp));
	    % Now need the longitudinal integrals with the normalization
	    % Most of them are zero in this degenerate case, which is
            % nice. NOW TAKING OUT THE LONGITUDINAL INTEGRAL IS ONLY GOOD
            % FOR THE AXISYMMETRIC CASE, BUT NEEDN'T EVEN ENTER THIS LOOP
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
	    	    
	    if xver==1
	      % Produce Maple string which should be zero
	      if m1>0
		str1=sprintf('sin(%i*phi)',m1);
	      else
		str1=sprintf('cos(%i*phi)',abs(m1));
	      end
	      if m2>0
		str2=sprintf('sin(%i*phi)',m2);
	      else
		str2=sprintf('cos(%i*phi)',abs(m2));
	      end
	      str3=sprintf('exp(%i*I*phi)',-mp);
	      % fprintf(fid,'%s\n',...
	      % sprintf('a:=a+evalf(abs(int(%s*%s*%s,phi=0..2*Pi))-%11.9f);',... 
	      %      str1,str2,str3,absI)...
	      %)
	    end
	    % Put in the sqrt(2) if needed
	    absI=absI*sqrt(2-[m1==0])*sqrt(2-[m2==0]);

	    % Now add in the colatitudinal integral if the above wasn't zero
	    if absI~=0
	      % This is the colatitudinal integral, in absolute value
	      % Watch the polynomial accuracy of the integration
	      % This makes a HUGE difference if you get it wrong
	      intgab=absI*abs(gausslegendre([-1 1],integrand,[2*L+lp]));
	    else
	      intgab=0;
	    end
	    
	    % Check this integral was done right
	    % Wait... if it's NOT zero, what does that mean? 
	    % It's a complicated set of conditions of course, reflecting
            % the triangle conditions etc. best exemplified by the other
            % methods of calculation
	    % Is there only one abs(mp) in this sequence for which this is
            % nonzero?
	     if intgab~=0
	     disp(sprintf(...
	     'a = %i b = %i lp = %i m1 = %+i m2 = %+i mp = %+i intgab = %8.3f',...
	     a,b,lp,m1,m2,mp,(2-[mp==0])*intgab))	      
	     end
	    
	    % Note that with lp=0 and mp=0 this must illustrate the orthogonality
	    % of the galm at fixed order
	    if xver==1 && mp==0 && lp==0
	      if a==b
		% The inner product of two identical functions, i.e. 1
		% Noting that ylm(0,0,0)=1/sqrt(4*pi)
		difer(intgab*sqrt(4*pi)-1)
	      else
		% The inner product of two orthogonal functions, i.e. 0
		difer(intgab) % Watch the indexing
	      end
	    end
	    % Don't forget to square
	    % Hold on: if you've shortcut mpcand to only be positive,
            % need to add this term twice, right?
	    % Cab(a,b)=Cab(a,b)+GC(index)*intgab^2;
            % This is if you never do any negatives and cut the loop in
            % almost half
	    Cab(a,b)=Cab(a,b)+(2-[mp==0])*GC(index)*intgab^2;
	  end
	end
      end
    end
    
    % Normalize, let the signal strength be S=1
    Cab=Cab/sqrt(pi);
   case 0 % SLOW, EVEN FOR THE TINIEST PROBLEMS
    Cab=repmat(0,Lpot,Lpot);
    for a=1:Lpot
      for b=a:Lpot
	% disp(sprintf('Tapers %i %i',a,b))
	for lp=0:2:min(2*l,2*L)
	  for mp=-lp:lp
	    intgab=0;
	    % This is a double sum already
	    for one=1:length(EMG)
	      l1=ELG(one);
	      m1=EMG(one);
	      %disp(sprintf('First loop %i %i',l1,m1))
	      % This is a double sum already
	      
	      % WHY WOULDN'T WE CUT THIS LOOP NOW?
	      % Attempt to cut loop, here, too, but leave uncut loop for
	      % verification; using or not using the selection rule
	      % noselrul=1:length(EMG);
	      % Put in all three of the selection rules
	      selrul=find(~mod(ELG+lp+l1,2) & ...
			  EMG==mp-m1 & ELG<=lp+l1 & ELG>=abs(lp-l1));
	      for two=selrul(:)'
		l2=ELG(two);
		m2=EMG(two);
		%disp(sprintf('Second loop %i %i',l2,m2))
		bprime=addmabout(L,ELGB(b),EMGB(b));
		% Any zeroes in here? Then you're doing too much work
		crap=GLMAB(one,a).*GLMAB(two,bprime)*...
		     sqrt(2*l1+1)*sqrt(2*l2+1)*...
		     indeks(wigner3jm(lp,l1,l2,0,0,0),'end')*...
		     indeks(wigner3jm(lp,l1,l2,-mp,m1,m2),'end');
		% If crap is zero, better be on account of GLMAB and not
		% WIGNER... maybe build in that we know ahead of time
                % which m we are dealing with.
		intgab=intgab+crap;
	      end
	    end
	    if xver==1 && lp==0
	      % This must illustrate the orthogonality; lp=0,
	      % sqrt(2*lp+1)=1; we didn't put in the sqrt(4*pi) so we don't
	      % have to compare it either; it is now either 0 or +/-1; it
	      % can be -1 since the sign is arbitrary for the eigenfunctions
	      if a==b
		difer(abs(intgab)-1)
	      else
		difer(intgab)
	      end
	    end
	    Cab(a,b)=Cab(a,b)+...
		     abs(intgab)^2*(2*lp+1)*...
		     indeks(wigner3jm(l,l,lp,0,0,0),'end')^2;
	  end
	end
      end
    end
    Cab=Cab/8/pi^2;
   otherwise
    error('Specify correct method')
  end

  % Symmetrize in one step
  Cab=Cab+tril(Cab',-1);

  % Better sort them in function of decreasing eigenvalue
  [V,j]=sort(V,'descend');

  % Since the bprime fix this is now properly sorted
  Cab=Cab(j,j);  

  % Is this then the global ordering, one would think so
  morder=EMGB(j);
  % You can get this from resorting GLMALPHA or straight from GALPHA 

  % Check the sum perpendicular to the diagonal, which is meaningless
  if xver==1
    for index=-Lpot+1:Lpot-1
      s(index+Lpot)=sum(diag(fliplr(Cab),index));
    end
  end
  eval(sprintf('save %s Cab l L TH meth morder',fnpl))
  disp(sprintf('Elapsed CPU time %8.3f',cputime-t))
end

% Put the correct normalization for comparison with GAMMAB which is
% correct
Cab=Cab*(4*pi)^2;


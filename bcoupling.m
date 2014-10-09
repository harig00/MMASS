function varargout=bcoupling(TH,Lmax,sord,meth,xver)
% [K,A]=BCOUPLING(TH,Lmax,sord,meth,xver)
%
% Calculates the coupling matrix due to a spherical boxcar.
% Dahlen & Simons (2008), Eq. (57)
%
% INPUT:
%
% TH       Colatitudinal radius of the CAP, in degrees
%          Colatitudinal halfwidth of the BELT, degrees
% Lmax     Maximum degree of this non-bandlimited function
% evens    0 For all degrees [default]
%          1 Only do the even degrees, e.g. for PERIODOVAR
% sord     1 Single cap of diameter 2TH [default]
%          2 Double cap left by subtracting belt of width 2TH
%          3 Equatorial belt of width 2TH
% meth     1 The standard way with the explicit Wigner 3j symbols [default]
%          2 Via the spatiospectral concentration matrix 
% xver     1 Provides extra verification
%          0 Doesn't [default]
%
% OUTPUT:
%
% K        The boxcar coupling matrix
% A        The area of the region in question
%
% EXAMPLE:
%
% bcoupling('demo') % Should return nothing
%
% Last modified by fjsimons-at-alum.mit.edu, 05/09/2007

if ~isstr(TH)

  defval('TH',30)
  defval('Lmax',20)
  defval('sord',1)
  defval('meth',1)
  defval('xver',0)

  % Needs change - should provide even if Lmax stored is too high
  fname=fullfile(getenv('IFILES'),'BCOUPLING',...
		 sprintf('BCOUPLING-%i-%i-%i-%i.mat',...
			 TH,Lmax,sord,meth));

  % Calculates the areas of these regions
  A=4*pi*spharea(TH,sord);

  % Initialize
  K=repmat(NaN,Lmax+1,Lmax+1);

  if exist(fname)==2
    disp(sprintf('load %s',fname))
    load(fname)
  else
    % Do the trivial case as well only if you are in verification mode
    if xver~=1
      if (TH==180 & sord==1) || (TH==0 & sord==2)
	K=eye(size(K));
	A=4*pi;
	save(fname,'K','Lmax','TH','sord','A','meth')
	% Output
	varns={K,A};
	varargout=varns(1:nargout);
	return
      end
    end
    switch meth
     case 1
      % Calculates the power spectrum of the boxcar up to the maximum
      % degree ever needed
      [Bl,dels]=bpboxcap(TH,2*Lmax,[],0,sord);
      
      % Get all the Wigner symbols at once up the the maximum degree ever
      % needed 
      [jk,C0,S0,LW]=zeroj(0,0,2*Lmax);

      h=waitbar(0);
      % Calculates the coupling kernel  
      for l=0:Lmax
	waitbar(l/Lmax,h,'BCOUPLING: Performing the calculations')
	for lp=l:Lmax
	  % Get the Wigner 3j symbols with bottom-row of zero
	  % W=wigner0j(Lmax,l,lp); % This is slow
	  % Get the Wigner 3j symbols with bottom-row of zero
	  W=zeroj(0:l+lp,l,lp,LW,[],C0,S0); % This is fast
	  
	  % Construct the symmetric part
	  K(l+1,lp+1)=sum((2*[0:l+lp]+1).*Bl(1:l+lp+1)'.*W.^2);
	  % Symmetrize
	  K(lp+1,l+1)=K(l+1,lp+1);
	end
      end
      close(h)
      % Don't forget about the (2l'+1) for the column dimensions
      K=K.*repmat(2*[0:Lmax]+1,Lmax+1,1)/A;
      % Or do K=K/A*diag(2*[0:Lmax]+1);
     case 2  
      % Calculate the localization kernel
      D=dlmlmp(TH,Lmax,Lmax,sord);
      % Compare with SDWCAP if you will - just check order 0
      if xver==1
	if sord==1
	  [E,V,N,th,C,ngl1,ngl2,unc,com,sdl,K]=sdwcap(TH,Lmax,0);
	elseif sord==2
	  [E,V,th,C,ngl1,ngl2,K]=sdwcap2(90-TH,Lmax,0);
	end
	difer(D(1:Lmax+1,1:Lmax+1)-K)
      end
      % Undo the block sorting and square
      [EM,EL,mz,blkm,dblk]=addmout(Lmax); 
      
      % Shouldn't I be making this complex now?
      %  U=ummp(Lmax); % This doesn't seem to help at all
      U=eye((Lmax+1)^2);
      % Square this thing pointwise
      D=(U*D(dblk,dblk)).^2;
      % Add all these orders together
      for l=0:Lmax
	b=l^2+1;
	e=(l+1)^2;
	% Construct the symmetric part
	for lp=l:Lmax
	  bp=lp^2+1;
	  ep=(lp+1)^2;
	  % The symmetric expression
	  K(l+1,lp+1)=sum(sum(D(b:e,bp:ep)));
	  K(lp+1,l+1)=K(l+1,lp+1);
	end
      end
      % The asymmetric expression
      K=K./repmat(2*[0:Lmax]'+1,1,Lmax+1)*4*pi/A;
     otherwise
      error('Specify valid method')
    end
    save(fname,'K','Lmax','TH','sord','A','meth')
  end

  % Check the row sum
  disp(sprintf('BCOUPLING Max row sum %6.2f %s',max(sum(K,2))*100,'%'))

  % When TH=180 should get the identity matrix I would think
  if (TH==180 & sord==1) || (TH==0 & sord==2)
    difer(K-eye(size(K)))
  end

  % Output
  varns={K,A};
  varargout=varns(1:nargout);

elseif strcmp(TH,'demo')
  % Note that Lmax is merely the size of the kernel; LBmax affects the
  % quality of the approximation. But LBmax has to be only 2*Lmax for there
  % to be no more difference upon increasing LBmax - by the selection
  % rules, regardless of where you are in the boxcar expansion...  even
  % though the kernel won't be complete and the rowsum not one... increases
  % in Lmax must go hand in hand with increases in LBmax. So LBmax is not
  % a free parameter here. This way there should be no difference in the
  % overlapping sections of two coupling kernels computed out to different
  % degrees. 
  [K2,A]=bcoupling(30,60,1,1); [K,A]=bcoupling(30,40,1,1);
  difer(K2(1:length(K),1:length(K))-K)
  [K1,A]=bcoupling(30,40,1,1); [K2,A]=bcoupling(30,40,1,2);
  difer(K2-K1)
  % Test the whole-sphere thing
  bcoupling(180,50)
end


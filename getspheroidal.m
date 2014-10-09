function varargout=getspheroidal(mod1,dir1,catflag)
% [rad,nn,el,ww,U,V,P,dUdr,dVdr,dPdr,qq]=GETSPHEROIDAL(mod1,dir1,catflag);
%
% Loads SPHEROIDAL normal-mode radial eigenfunctions and frequencies.
%
% INPUT:
%
% mod1        The model file name [default: 'prem800']
% dir1        The directory name [default: $IFILES/EARTHMODELS/MODES/blabla]
% catflag     0 Modes files made by YANNOS aka MINOSY [default]
%             1 Modes files made by MINOS_BRAN and formatted by
%               EIGCON from the packages MINEOS-1.0.0 or MINEOS-1.0.2; 
%               the former is a small catalog, the latter a big
%               one... but remember the Earth models are slightly different!
% 
% OUTPUT:
%
% rad         Radii at which the radial eigenfunctions are evaluated [m]
% nn          Array of branch numbers
% el          Array of spherical harmonic degrees (sorted)
% ww          Array of angular frequencies [rad/s]
% U, dUr      Radial radial eigenfunctions and derivatives
%             Although "displacement", units are [kg^{-1/2}, kg^{-1/2}/m]
% V, dVr      Tangential radial eigenfunctions and derivatives 
%             Although "displacement", units are [kg^{-1/2}, kg^{-1/2}/m]
% P, dPdr     Potential-perturbation radial eigenfunctions and derivatives
%             Although "potential", units are [m/s^2/sqrt(kg) 1/s^2/sqrt(kg)]
% qq          Modal decay rate, "gamma" defined by DT (9.53)
%
% EXAMPLES:
%
% getspheroidal('demoX',[],cf) where X=1 through 9 and cf=0 or 1
% These should reproduce figures from Dahlen & Tromp (DT) but with gravity
% potential eigenfunctions included also. Note that between catalogs, and
% in comparison to the DT figures, certain modes may have been
% mislabeled, signs may be reversed, etc.
% 
% SEE ALSO: GETTOROIDAL, MODESTRAIN, EQPOTENTIAL
%
% NOTES:
%
% For notation, see Dahlen & Tromp (DT), p. 289. But be aware that
% their functions V are a factor sqrt(l(l+1)) bigger than ours.
% 
% Last modified by jchawtho-at-princeton.edu, 07/01/2007
% Last modified by fjsimons-at-alum.mit.edu, 07/16/2010
% Last modified by efwelch-at-princeton.edu, 07/11/2012

% Gravitational constant
G=fralmanac('GravCst','Earth');

defval('catflag',0)

defval('mod1',[])

if isempty(strmatch('demo',mod1))
  % Extra normalization verification
  defval('xver',0)
  
  % Progress string
  stronk='Read %6i modes in %2i s at an average of %5.3f s per kmode';

  switch catflag
   
   case 0
    % YANNOS produced catalogs

    % Get directory and model name, and construct modes file name 
    defval('mod1','prem800')
    defval('dir1',fullfile(getenv('IFILES'),'EARTHMODELS','MODES','YANNOS'))
    fname=fullfile(dir1,sprintf('fct_%sS.direct',mod1));
    disp(sprintf('Loading %s',fname))
    
    % Prepare for output in a format understandable to human intelligence
    [radmod,rho]=getmodel(mod1,dir1,catflag);
    
    % Get the file size to be smart about initialization
    fs=fsize(fname);
    
    % %%%%%%%%%%
    % Turns out this takes a whole long while to load
    fid=fopen(fname);
    % This is the length of the header and all other records
    len=fread(fid,1,'int32');
    % Why -1, no idea
    minus1=fread(fid,2,'int32');
    % Number of ocean layers in the model
    noc=fread(fid,1,'int32');
    % Number of layers used
    nrad=fread(fid,1,'int32');
    % Radii divided by the Earty radius
    % Note that the radii and the eigenfunctions go in different directions
    rd=fread(fid,nrad,'float32');
    % Presumably all zero - so the header has the same format as the data
    trash=fread(fid,len/4-nrad-5,'int32');
    % Check this is indeed the case
    difer(trash,[],[],NaN)
    
    % The part above here is technically the file header
    difer(len-ftell(fid),[],[],NaN)
    
    % What is, ultimately, the number of modes read in?
    nmodes=fs/len-1;
    
    disp(sprintf('Preparing to read %i SPHEROIDAL modes',nmodes))
    
    % Initializing arrays saves a LOT of time here
    [nn,el,ww,qq,cg]=deal(zeros(nmodes,1));
    d=zeros(nrad,6,nmodes);
    
    % On to the "read" of the modes file proper
    nxt=fread(fid,1,'int32');
    i=0; tic
    while ~isempty(nxt)
      % With every increment you're moving by 'len' positions in the file
      i=i+1;

      % Overtone "branch" number
      nn(i)=nxt;

      % Read the mode-identifying header bits
      [el(i),ww(i),qq(i),cg(i)]=readbitsy(fid);
      
      % Eigenfunctions and radial derivative: pU,dUdr,V,dVdr,P,dPdr
      d(1:nrad,1:4,i)=[reshape(fread(fid,4*nrad,'float32'),[4 nrad])]';
      % Could possibly combine even these two lines here
      d(1:nrad,5:6,i)=reshape(fread(fid,nrad*2,'float32'),nrad,2);
       
      % Go on and report on progress
      nxt=moveon(fid,i,stronk);
    end
    
    % Done
    disp(sprintf(stronk,i,round(toc),toc/i*1000)); 
    fclose(fid);

    % These are ALMOST but not quite identical to the quoted values in mod1 
    radius=radmod(end);
    rad=rd*radius;
    
    % Report on this annoying misfit
    difer([rad-radmod]/radius,5,[],NaN)
    
    % Figure out the normalization factor
    [norma,rhomean]=rhobar(rad,rd,rho);
    
    % Sort accordingly and normalize
    d=flipdim(d,1)/sqrt(norma);
    
    % Make it right
     U  =squeeze(d(:,1,:));
    dUdr=squeeze(d(:,2,:))/radius;
     V  =squeeze(d(:,3,:));
    dVdr=squeeze(d(:,4,:))/radius;
    
    % P is normalized slightly differently
     P  =squeeze(d(:,5,:)).*pi.*rhomean.*G.*radius;
    dPdr=squeeze(d(:,6,:)).*pi.*rhomean.*G;
    % Note this is close to 1/rhomean*10000*4*pi/3
    
    % This changes qq into the exponential decay term "gamma" from DT (9.53)]
    % which is more useful and is also easily extractable from MINEOS catalogs.
    qq=ww./qq/2;
   
    clear d
   case 1
    % MINEOS produced modes (after formatting by EIGCON)
    
    % Get directory and model name, and construct modes file names
    % An anisotropic prem
    defval('mod1','aniprem808')
    defval('dir1',fullfile(getenv('IFILES'),'EARTHMODELS','MODES','MINEOS-1.0.0'))
    % An anisotropic prem with the ocean removed
    defval('mod1','noocean_aniprem808')
    defval('dir1',fullfile(getenv('IFILES'),'EARTHMODELS','MODES','MINEOS-1.0.2'))
    % Note that for this latter case, the "tightened" eigenfunction files
    % cannot be read by this program. Evan moved beyond this.     
    
    % Radial modes
    fnameR=fullfile(dir1,sprintf('%s_R',mod1));
    % Spheroidal modes
    fnameS=fullfile(dir1,sprintf('%s_S',mod1));
    disp(sprintf('Loading %s',fnameR))
    disp(sprintf('Loading %s',fnameS))

    % Prepare for output in a format understandable to human intelligence
    [radmod,rho]=getmodel(mod1,dir1,catflag);
    nrad=length(radmod);
    
    % Get the file size (bytes) to be smart about initialization
    fsS=fsize(fnameS);
    fsR=fsize(fnameR);
    
    % Compute number of modes stored knowing each word is 4 bytes
    nmodesS=fsS/28/(1+nrad);
    nmodesR=fsR/4/(7+3*nrad);
    nmodes=nmodesR+nmodesS;
    
    % Initializing arrays saves a LOT of time here
    [nn,el,ww,qq,rn,vn,an]=deal(zeros(nmodes,1));
    d=zeros(nrad,7,nmodes);

    % Open binary files
    fidS=fopen(fnameS);
    fidR=fopen(fnameR);
    disp(sprintf('Preparing to read %i SPHEROIDAL + RADIAL modes',nmodes))

    % On to the "read" of the RADIAL modes file itself
    i=0; tic;
    nxt=fread(fidR,1,'int32');
    while ~isempty(nxt)
      % With every increment you're moving by (nrad*3+7) positions in the file
      i=i+1;

      % Overtone "branch" number
      nn(i)=nxt;

      % Read the mode-identifying header bits
      [el(i),ww(i),qq(i),rn(i),vn(i),an(i)]=readbitsm(fidR);

      % Read radii and eigenfunctions
      d(:,1:3,i)=reshape(fread(fidR,nrad*3,'float32'),3,nrad)';

      % Go on and report on progress
      nxt=moveon(fidR,i,stronk);
    end

    % On to the "read" of the SPHEROIDAL-proper modes file itself
    nxt=fread(fidS,1,'int32');
    while ~isempty(nxt)
      % With every increment you're moving by (nrad*7+7) positions in the file
      i=i+1;

      % Overtone "branch" number
      nn(i)=nxt;
      
      % Read the mode-identifying header bits
      [el(i),ww(i),qq(i),rn(i),vn(i),an(i)]=readbitsm(fidS);
      
      % Read radii and eigenfunctions
      d(:,:,i)=reshape(fread(fidS,nrad*7,'float32'),7,nrad)';

      % Go on and report on progress
      nxt=moveon(fidS,i,stronk);
    end

    % Done
    disp(sprintf(stronk,i,round(toc),toc/i*1000)); 
    fclose(fidS); fclose(fidR);

    % Check that the prediction on the number of modes was right
    difer(length(el)-nmodes,[],[],NaN)
    
    % Double check that these fields are in fact homogeneous
    difer(sum(rn~=rn(1) | vn~=vn(1) | an~=an(1)),[],[],NaN)

    % Now every normalization is unique
    rn=rn(1); vn=vn(1); an=an(1);

    % Sort by degree
    [el,indices]=sort(el);
    nn=nn(indices);
    ww=ww(indices);
    qq=qq(indices);
    d=d(:,:,indices);
    
    % Keep just one of the modes radii
    rd=squeeze(d(:,1,1));
    
    % These are ALMOST but not quite identical to the quoted values in mod1 
    radius=max(radmod);
    rad=rd*radius;
    
    % Report on this annoying misfit which under single precision is
    % worse than for YANNOS
    difer([rad-radmod]/radius,1,[],NaN)
    
    % Figure out the normalization factor
    [norma,rhomean]=rhobar(rad,rd,rho);

    % Now normalize to D&T convention (realizing that their V and W's
    % are sqrt(l*(l+1)) bigger than ours
    d=d/sqrt(norma)/vn*rn;

    % Include the dependence of eigenfunctions on frequency ww
    bigW=repmat(ww',nrad,1);
    % Make it right
     U  =squeeze(d(:,2,:)).*bigW;
    dUdr=squeeze(d(:,3,:)).*bigW/radius;
     V  =squeeze(d(:,4,:)).*bigW;
    dVdr=squeeze(d(:,5,:)).*bigW/radius;
    
    % P is normalized slightly differently
     P  =squeeze(d(:,6,:)).*bigW.*pi.*rhomean.*G.*radius;
    dPdr=squeeze(d(:,7,:)).*bigW.*pi.*rhomean.*G;
    
    clear d bigW
  end
  
  % Verify that these are "normal" modes as in DT (8.107)
  if xver==1
    % Number of tests to run at random
    inmax=25;    

    % initialize to save time (these are used in the P verifications)
    [A,B]=deal(zeros(length(rad),1));    
    % The current MINEOS-1.0.0 catalog has an abundance of ignored P
    % modes due to the way the minos_bran paramter, wgrav was
    % set. Essentially P modes with frequencies above 10 mHz are
    % not computed. The loop below will only normalize modes that
    % are nonzero. This could actually be applied to both catalogs
    % for good measure. 
    these=find(sum(P,1)~=0);
    nonzeros=length(these);
    disp(sprintf('This catalog has %i non zero P modes',nonzeros))

    [ppsave,nplsave]=deal(zeros(length(rad),inmax));
    for in=1:inmax
      % Test and report on the inner product normalization
      somod=ceil(rand(2,1)*size(U,2));
      innerproducts(U,V,nn,el,rad,rho,somod);

      % Test only "these" that are gravitationally affected
      somod=these(ceil(rand(1,1)*length(these)));
      
      % Rename the particular mode to improve clarity
      zel=el(somod); zen=nn(somod);
      nUl=U(:,somod); nVl=V(:,somod); nPl=P(:,somod);
      
      % Now begin check of P using DT eq. 8.55 at every radius
      for jind=1:length(rad)
	range1=1:jind;
	range2=jind:length(rad);
	rr=rad(jind);
	% First term in (8.55) 
	if jind>1
	  A(jind)=trapeze(rad(range1),(rad(range1)./rr).^(zel+1).*rho(range1).* ...
			  zel.*(nUl(range1)+(zel+1).*nVl(range1)));
	end
	% Second term in (8.55) 
	if jind~=length(rad)
	  B(jind)=trapeze(rad(range2),rho(range2).*(rr./ ...
	    rad(range2)).^(zel).*(-(zel+1).*nUl(range2)+zel*(zel+1).* ...
				  nVl(range2)));
	end
      end
      
      % Save what the P's should be from DT 8.55 
      ppsave(:,in)=-((4*pi*G)/(2*zel+1))*(A+B);

      % Save what they are from YANNOS/MINEOS
      nplsave(:,in)=nPl;

      % Save which mode we used
      whichmod(:,in)=[zen; zel];
    end
   
    % Don't consider NaN's that arose from integration
    % Should be array of ones
    response=input(['Plot comparison between expected P and actual', ...
		    ' P for random modes? (yes or no)'],'s');
    if strcmp(response,'yes')
      for jj=1:inmax
	clf
	thismod=ceil(size(ppsave,2)*rand(1,1));
	plot(ppsave(:,thismod),rad,'r');
	hold on; 
	plot(nplsave(:,thismod),rad);
	legend('Eq. 8.55',sprintf('P, n=%d, l=%d',whichmod(1,thismod),...
				  whichmod(2,thismod)));
	pause
      end
    end
    % Check this to see there is no bias
    % plot(ppsave(:)./nplsave(:),'.')
    
    % Put in some testing for qq possibly?  E.g. (9.54) 
  end
  
  % Check the guess of the number of eigenvectors was good
  difer(nmodes-size(U,2),[],[],NaN)
  
  % Check the degrees were sorted
  difer(el-sort(el),[],[],NaN)
  
  % Prepare output
  vars={rad,nn,el,ww,U,V,P,dUdr,dVdr,dPdr,qq};
  varargout=vars(1:nargout);  
else
  % It is one of the demos!
  demnum=mod1;
  [rad,nn,el,ww,U,V,P]=getspheroidal([],[],catflag);
  radius=rad(end);
  if strcmp(demnum,'demo1')
    % Should be close to DT Figure 8.14
    % ScSv modes
    ens=[9 12 14 21 23 26 30 33 35 42];
    els=repmat(2,size(ens));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo2')
    % Should be close to DT Figure 8.14
    % PKIKP modes
    ens=[8 10 13 18 22 25 27 29 31 34];
    els=repmat(2,size(ens));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo3')
    % Should be close to DT Figure 8.14
    % JSv modes (Note the likely typo)
    ens=[6 11 15 16 19 24 28 32 40 45];
    els=repmat(2,size(ens));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo4')
    % Should be close to DT Figure 8.15
    % CMB Stoneley modes (Note the likely typo)
    ens=[1 2 3 4 5 6 6 7];
    els=10:10:90;
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo5')
    % Should be close to DT Figure 8.15
    % ICB Stoneley modes
    ens=[2 3 4 5 6 7 8 9];
    els=[2 4 6 9 11 14 18 21];
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo6')
    % Should be close to DT Figure 8.11
    % Fundamental Rayleigh modes
    els=30:10:130;
    ens=zeros(size(els));
    ylimo=[radius-1200000 radius];
  elseif strcmp(demnum,'demo7')
    % Should be close to DT Figure 8.11
    % First overtone Rayleigh modes
    els=30:10:130;
    ens=ones(size(els));
    ylimo=[radius-1200000 radius];
  elseif strcmp(demnum,'demo8')
    % Should be close to DT Figure 8.11
    % Second overtone Rayleigh modes
    ylimo=[radius-1200000 radius];
    els=30:10:130;
    ens=repmat(2,size(els));
  elseif strcmp(demnum,'demo9')
    % Should be close to DT Figure 8.16
    % First ten radial modes
    ylimo=[radius-6.371e6 radius-600e3];
    ens=0:9;
    els=zeros(size(ens));
  else
    error('Specify valid demo string')
  end

  % Now move to making the actual plot
  ylabs=sort(radius-[670 2899 5149.5]*1000);
  clf
  [ah,ha,H]=krijetem(subnum(1,length(ens)));
  for inx=1:length(ens)
    axes(ah(inx))
    witsj=nn==ens(inx) & el==els(inx);
    up(inx)=plot(U(:,witsj),rad,'-'); hold on
    pp(inx)=plot(P(:,witsj),rad,'k-'); 
    vp(inx)=plot(sqrt(els(inx)*(els(inx)+1))*V(:,witsj),rad,'--');
    axis tight
    plot(xlim,repmat(ylabs,2,1),':')
    hold off
    t(inx)=title(sprintf('_{%i}S_{%i}',ens(inx),els(inx)));
    ylim(ylimo)
    drawnow
  end
  % Cosmetics
  set(ah,'ytick',ylabs,...
	 'xtick',[],'ytickl',{'ICB','CMB',670})
  set(ah(2:end),'visible','off')
  set(cell2mat(gettit(ah(2:end)))','visible','on')
  set(ah(1),'color','none','box','off')
  longticks(ah,1)
  nolabels(ah(2:end),2)
  nolabels(ah,1) 
  shrink(ah,1/1.5,1)    
  fig2print(gcf,'portrait')
  figdisp([],sprintf('%s_%i',demnum,catflag),[],1)
  % If it's ugly, run it again!
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [el,ww,qq,cg]=readbitsy(fid)
% Spherical harmonic degree
el=fread(fid,1,'int32');
% Eigenfrequency
% Units are Hz times 2*pi
ww=fread(fid,1,'double');
% Modal quality factor (DT 9.54)
qq=fread(fid,1,'double');
% Group velocity FJS???
cg=fread(fid,1,'double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [el,ww,qq,rn,vn,an]=readbitsm(fid)
% Spherical harmonic degree
el=fread(fid,1,'int32');
% Eigenfrequency
% Units are Hz times 2*pi
ww=fread(fid,1,'float32');
% Modal gamma (DT 9.53)
qq=fread(fid,1,'float32');
% Radius normalization,
rn=fread(fid,1,'float32');
% Velocity normalization,
vn=fread(fid,1,'float32');
% Acceleration normalization
an=fread(fid,1,'float32');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nxt=moveon(fid,i,stronk)
nxt=fread(fid,1,'int32');
% Advise on progress
if ~mod(i,3000)
  t=toc; disp(sprintf(stronk,i,round(t),t/i*1000)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function innerproducts(U,V,nn,el,rad,rho,somod)
% String to report
forms='Normalization at n = %3i; np = %3i ; l = %3i; lp = %3i';
% The relevant degrees and branch numbers
zel =el(somod(1));  zen=nn(somod(1));
zelp=el(somod(2)); zenp=nn(somod(2));
% The relevant radial eigenfunctions
nUl=U(:,somod(1)); nUlp=U(:,somod(2));
nVl=V(:,somod(1)); nVlp=V(:,somod(2));
% The inner products that define the normalization
sb1=trapeze(rad,rad.^2.*rho.*(nUl.^2+zel.*(zel+1).*nVl.^2));
sb2=trapeze(rad,rad.^2.*rho.*(nUlp.^2+zelp.*(zelp+1).*nVlp.^2));
sb0=trapeze(rad,rad.^2.*rho.*(nUl.*nUlp+...
			      sqrt(zel.*(zel+1)).*...
			      sqrt(zelp*(zelp+1)).*nVl.*nVlp));
% Report on the results
disp(sprintf(...
    [forms ' %8.3f %8.3f %8.3f'],...
    zen,zenp,zel,zelp,abs(sb1),abs(sb2),abs(sb0)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [norma,rhomean]=rhobar(rad,rd,rho);
% The mean density of the Earth
rhomean=3*trapeze(rd,rd.^2.*rho);
% This is the density of the Earth when all the mass is squeezed into a
% sphere of unit radius
norma=3*trapeze(rad,rad.^2.*rho);

disp(sprintf('Mean model density: %8.3f',rhomean))

% Note that these were hardwired in both codes
rhomean=5515; 
norma=max(rad)^3*5515;

disp(sprintf('Mean model density used: %8.3f',rhomean))


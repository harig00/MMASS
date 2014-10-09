function varargout=gettoroidal(mod1,dir1,catflag)
% [rad,nn,el,ww,W,dWdr,qq]=GETTOROIDAL(mod1,dir1,catflag);
%
% Loads TOROIDAL normal-mode radial eigenfunctions and frequencies.
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
% W, dWdr     Toroidal radial eigenfunctions and derivatives
%             Although "displacement", units are [kg^{-1/2}, kg^{-1/2}/m]
% qq          Modal decay rate, "gamma" defined by DT (9.53)
%
% EXAMPLES:
%
% gettoroidal('demoX',[],cf) where X=1 through 6 and cf=0 or 1
% These should reproduce figures from Dahlen & Tromp (DT). Note that
% between catalogs, and in comparison to the DT figures, certain modes
% may have been mislabeled, signs may be reversed, etc.
% 
% SEE ALSO: GETSPHEROIDAL, MODESTRAIN, EQPOTENTIAL
%
% NOTES:
%
% For notation, see Dahlen & Tromp (DT), p. 289. But be aware that
% their functions W are a factor sqrt(l(l+1)) bigger than ours.
% 
% Last modified by jchawtho-at-princeton.edu, 07/01/2007
% Last modified by fjsimons-at-alum.mit.edu, 07/16/2010
% Last modified by efwelch-at-princeton.edu, 07/11/2012

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
    fname=fullfile(dir1,sprintf('fct_%sT.direct',mod1));
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
    [nn,el,ww,qq,cg]=deal(nan(nmodes,1));
    d=nan(nrad,2,nmodes);
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
      
      % Eigenfunction and radial derivative: W and dWdr
      d(1:nrad,1:2,i)=[reshape(fread(fid,2*nrad,'float32'),[2 nrad])]';
       
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
     W  =squeeze(d(:,1,:));
    dWdr=squeeze(d(:,2,:))/radius;
    
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

    fname=fullfile(dir1,sprintf('%s_T',mod1)); % Toroidal modes
    disp(sprintf('Loading %s',fname))

    % Prepare for output in a format understandable to human intelligence
    [radmod,rho]=getmodel(mod1,dir1,catflag);
    nrad=length(radmod);
    
    % Get the file size (bytes) to be smart about initialization
    fs=fsize(fname);

    % Compute number of modes stored knowing each word is 4 bytes
    nmodesT=fs/4/(7+3*nrad);
    nmodes=nmodesT;
    
    % Initializing arrays saves a LOT of time here
    [nn,el,ww,qq,rn,vn,an]=deal(nan(nmodesT,1));
    d=nan(nrad,3,nmodesT);

    % Open binary files
    fid=fopen(fname);
    disp(sprintf('Preparing to read %i TOROIDAL modes',nmodesT))

    % On to the "read" of the TOROIDAL modes file itself
    i=0; tic;
    nxt=fread(fid,1,'int32');
    while ~isempty(nxt)
      % With every increment you're moving by (nrad*3+7) positions in the file
      i=i+1;

      % Overtone "branch" number
      nn(i)=nxt;

      % Read the mode-identifying header bits
      [el(i),ww(i),qq(i),rn(i),vn(i),an(i)]=readbitsm(fid);
      
      % Read radii and eigenfunctions
      d(:,:,i)=reshape(fread(fid,nrad*3,'float32'),3,nrad)';

      % Go on and report on progress
      nxt=moveon(fid,i,stronk);
    end

    % Done
    disp(sprintf(stronk,i,round(toc),toc/i*1000)); 
    fclose(fid);

    % Check that the prediction on the number of modes was right
    difer(length(el)-nmodesT,[],[],NaN)
    
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
     W  =squeeze(d(:,2,:)).*bigW;
    dWdr=squeeze(d(:,3,:)).*bigW/radius;
    
    clear d bigW
  end

  % Verify that these are "normal" modes as in DT (8.107)
  if xver==1
    % Number of tests to run at random
    inmax=25;    
    
    for in=1:inmax
      % Test and report on the inner product normalization
      somod=ceil(rand(2,1)*size(W,2));
      innerproducts(W,nn,el,rad,rho,somod)
    end
  end

  % Check the guess of the number of eigenvectors was good
  difer(nmodes-size(W,2),[],[],NaN)
  
  % Check the degrees were sorted
  difer(el-sort(el),[],[],NaN)
  
  % Prepare output
  vars={rad,nn,el,ww,W,dWdr,qq};
  varargout=vars(1:nargout);  
else
  % It is one of the demos!
  demnum=mod1;
  [rad,nn,el,ww,W,dWdr]=gettoroidal([],[],catflag);
  radius=rad(end);
  if strcmp(demnum,'demo1')
    % Should be close to DT figure 8.4
    % First ten fundamental modes
    els=2:11;
    ens=zeros(size(els));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo2')
    % Should be close to DT figure 8.6
    % First ten degree two modes
    ens=0:9;
    els=repmat(2,size(ens));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo3')
    % Should be close to DT figure 8.7
    % Fundamental Love wave equivalent modes
    els=[30 40 50 60 70 80 90 100 110 120 130];
    ens=repmat(0,size(els));      
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo4')
    % Should be close to DT figure 8.7
    % First-overtone Love wave equivalent modes
    els=30:10:130;
    ens=ones(size(els));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo5')
    % Should be close to DT figure 8.7
    % Second-overtone Love wave equivalent modes
    els=30:10:130;
    ens=repmat(2,size(els));
    ylimo=[0 radius];
  elseif strcmp(demnum,'demo6')        
    % Should be close to DT figure 8.7
    % Second-overtone Love wave equivalent modes
    ens=[13 12 11 10 9 8 7 6 5 4];
    els=[7 25 34 40 43 47 51 55 61 67];   
    ylimo=[0 radius];
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
    wp(inx)=plot(W(:,witsj),rad,'-'); hold on
    axis tight
    plot(xlim,repmat(ylabs,2,1),':')
    hold off
    t(inx)=title(sprintf('_{%i}T_{%i}',ens(inx),els(inx)));
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
function innerproducts(W,nn,el,rad,rho,somod)
% String to report
forms='Normalization at n = %3i; np = %3i ; l = %3i; lp = %3i';
% The relevant degrees and branch numbers
zel =el(somod(1));  zen=nn(somod(1));
zelp=el(somod(2)); zenp=nn(somod(2));
% The relevant radial eigenfunctions
nWl=W(:,somod(1)); nWlp=W(:,somod(2));
% The inner products that define the normalization
sb1=trapeze(rad,rad.^2.*rho.*(zel.*(zel+1).*nWl.^2));
sb2=trapeze(rad,rad.^2.*rho.*(zelp.*(zelp+1).*nWlp.^2));
sb0=trapeze(rad,rad.^2.*rho.*(sqrt(zel.*(zel+1)).*...
			      sqrt(zelp*(zelp+1)).*nWl.*nWlp));
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

% Note that thes were hardwired in both codes
rhomean=5515; 
norma=max(rad)^3*5515;

disp(sprintf('Mean model density used: %8.3f',rhomean))


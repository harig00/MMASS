function varargout=loris5(wav,N,J,precon,depkm,xver,gnorjr,panel)
% [vw,depkm,vstats]=LORIS5(wav,N,J,precon,depkm,xver,gnorjr,panel)
%
% Study of the sparsity in the models of Montelli and Ritsema
% Learn from the depth-dependence of scales in a seismic model.
%
% INPUT:
%
% wav      The chosen wavelet construction used, 'D2', 'D4' or 'D6'
% N        The power of the dyadic subdivision for the cubed sphere
% J        Maximum scale in both directions for the wavelet analysis
% precon   [1|0 1|0] toggles preconditioning in either dimension
% depkm    The single depth in km where the model is read and the results of the
%          calculations are saved - must match exactly. If a vector or
%          empty then you work on the result and get a plot
% xver     For extra verification showing some intermediate results
% gnorjr   1 For the model by Montelli et al. (2006)
%          2 For the model by Ritsema et al. (2010)
% panel    0 Does all panels [default]
%          1-6 Picks a specific panel
%
% OUTPUT:
%
% vw       The wavelet-transformed single-depth cubed-sphere for all chunks
%          or: the panel number if you're in all-depths output mode
% depkm    The depth in km, if you didn't already know it 
% vstats   A cell structure with the statistics per scale and depth
%
% EXAMPLE:
%
% panel=0;
% loris5('D4',7,4,[1 1],[],0,1,panel) % Paper figure A
% loris5('D4',7,4,[1 1],[],0,2) % Paper figure B
%
% Last modified by fjsimons-at-alum.mit.edu, 2/28/2011

% Identify defaults
defval('wav','D4')
defval('N',7)
defval('J',N-3-strcmp(wav,'D6'))
defval('precon',[1 1]);
% d=readGNmodel(NaN);
% This should hold both for Montelli and Ritsema
d=load('/u/fjsimons/IFILES/EARTHMODELS/MONTELLI/GNdepths');
defval('depkm',d/1000);
defval('xver',0);
defval('fs',11);
defval('panel',0);
% If 1, GN model, if 2, JR model
defval('gnorjr',1)

% Hardwired Earth radius, not ideal, see readGNmodel
rearthkm=6371;

% How many depth layers?
ndeps=length(depkm);

% The default is to do all depths and save the transform calculations if
% they don't exist for the next time you run the program
if gnorjr==1
  % And now define what will come out of the analysis at the end
  fname2=fullfile(getenv('IFILES'),'EARTHMODELS','MONTELLI','LORIS',...
		 sprintf('loris5_GN_%i_%i_%s_%i_%i_%i.mat',...
			 N,J,wav,precon(1),precon(2),panel));
elseif gnorjr==2
  % And now define what will come out of the analysis at the end
  fname2=fullfile(getenv('IFILES'),'EARTHMODELS','RITSEMA','LORIS',...
		 sprintf('loris5_JR_%i_%i_%s_%i_%i_%i.mat',...
			 N,J,wav,precon(1),precon(2),panel));
end

% If all of the depths are to be revisited
if ndeps>1
  if exist(fname2,'file')==2
    load(fname2)
  else 
    % The percentiles that will be CALCULATED
    percs=[85 100];
    % Study the properties of this one thresholding percentage
    tperc=85;
    % Title for the x axis
    xlstring1=sprintf('%s wavelet | scaling coefficients',wav);
    
    % Identify the scales in the mr=1 wavelet transform
    % Do one more than J so the last one is the scaling coefficients 
    [vwlev,vwlevs]=cube2scale(N,[J J]+1,1);
    
    % Potential check at a random depth
    deptsj=max(1,randi(ndeps));
    
    % Initialize the structure array for good measure
    fields={'deptkm','percts','mjuvar','totnum','tcount','normsq'};
    % And how many here will populate this?
    nmflds={1,    length(percs),2,    1,  1,   1   };
    nmcols=[1     J+1          J+1   J+1 J+1  J+1  ];
    nmrows=[ndeps ndeps       ndeps   1 ndeps ndeps];
    % Note that being this specific saves memory but increases speed
    vwstats=structempty(fields);
    for index=1:length(fields)
      vwstats.(fields{index})=...
	  cellnan([nmrows(index) nmcols(index)],1,nmflds{index});
    end
    % End initialization of the structure array with the results
    
    % Then we call ourselves once per depth layer, still for the straight J
    for index=1:ndeps
      % Remember that this function calls itself
      [vw,depchk,vstats]=...
	  loris5(wav,N,J,precon,depkm(index),xver,gnorjr,panel);
      difer(depchk-depkm(index),[],[],NaN)
      
      % Here we do the panel nulling, or should we do picking
      if ~~panel
	vw(:,:,skip(1:6,panel))=0;
	% Let's see what picking does instead
	vw=vw(:,:,panel);
	% Remember the vwlevs was just copies in 3D
	vwlevs=vwlev;
      end
      
      % Quick plot?
      if xver==1 && index==deptsj
	for jndex=1:J+1
	  vwnul=vw;
	  vwnul(vwlevs==jndex)=100;
	  title(sprintf('Check at depth %i for scale %i',index,jndex))
          % This subfunction makes a quick wavelet coefficient plot
	  clf
	  [cb,xcb]=quickplot(N,J,vwnul,depkm(index));
	  pause(1)
	end
      end
      
      % Remember, we work with the absolute values! But save the originals
      vworg=vw;
      vw=abs(vw);
      
      % Prepare for the study of the effects of their thresholding which is
      % here carried out per individual depth 
      kilme=vw<prctile(vw(:),tperc);
      nkept=sum(~kilme(:));
      
      % Take a look at what is being killed
      if xver==1
	for lndex=1:min(6,size(vw,3))
	  clf
	  spy(kilme(:,:,lndex)); hold on; fridplotw(N,J); axis off
	  title(num2str(lndex))
	  pause(1)
	end
      end
      
      % The overall norm upon inverse transformation is recovered as can be
      % checked, but we'll go straight on to identifying the
      % scale-dependent information and its inverse transformation
      normscsq=0;
      for jndex=1:J+1
	% Identify where the scales are located over ALL chunks (vwlevs!)
	itshere=[vwlevs==jndex];
	vwthere=vw(itshere);
	
	% Perform the reconstruction of the model using just this scale
	vwnul=vworg;
	vwnul(~itshere)=0;
	
	% Now calculate the inverse wavelet transforms, still for straight J
	switch wav
	 case 'D2'
	  % The 2-tap Haar/Daubechies wavelet transform
	  vwrec=angularD2WT(vwnul,[J J],'inverse',1);
	 case 'D4'
	  % The 4-tap Daubechies wavelet transform
	  vwrec=angularD4WT(vwnul,[J J],precon,'inverse',1);
	 case 'D6'
	  % The 6-tap Daubechies wavelet transform
	  vwrec=angularD6WT(vwnul,[J J],precon,'inverse',1);
	end
	
	% Let us not forget that the norm of the preconditioned transforms
	% changes, but when they are subjected to the inverse, the total
	% norm aligns itself with the original data, i.e. vstats.normv.
	% However, if there was no preconditioning, we can do better:
	vwrecnorm=norm(vwrec(:));
	if prod(precon)==0
	  difer(vwrecnorm-norm(vwnul(:)),[],[],NaN)
	  normscsq=normscsq+vwrecnorm^2;
	end
	
	% What is the scale-dependent percentage of the model mean square
	% signal strength that is being explained by the wavelet
	% coefficients at that scale? Note that the sum of the rms never
	% adds up to one hundred percent, and that the mean squared values
	% per scale only add to 100 if the transform is NOT preconditioned,
	% which provides the most insight. So we deal with the norm squared. 
	vwstats.normsq{index,jndex}=100*vwrecnorm^2/vstats.normv^2;
	
	% Quick plot?
	if xver==1
	  dax=[];
	  clf; ah=gca;
	  if ~~panel
	    % If it was only one panel, refake 3D
	    fak3=zeros([size(vwrec) 6]);
	    fak3(:,:,panel)=vwrec;
	    vwrec=fak3;
	  end
	  [cb,xcb,pgw]=quickplot(N,J,vwrec,depkm(index),dax);
	  for ons=1:6; delete(cat(1,pgw{ons}{:})); end
	  axes(ah); plotcont([],[],9)
	  % Also try this and make sure that it looks nice:
	  % plotoncube(vwrec,'2D',1,[],[],[],dax,[],0); 
	  % plotcont([],[],9)
	  % Also try this and make sure that it looks nice:
	  % plotoncube(vwrec,'3D',1,[],[],[],dax,[],0);
	  % shading flat; axis image; hold on; plotcont([],[],3)
	end
	
	% Make the diagnostics vector, row depth column scale
	vwstats.percts{index,jndex}=prctile(vwthere,percs);
	vwstats.mjuvar{index,jndex}=[mean(vwthere) var(vwthere)];
	
	% Keep the information in the array for good measure
	if jndex==1
	  vwstats.deptkm{index}      =depkm(index);
	end
	if index==1
	  vwstats.totnum{      jndex}=sum(itshere(:));
	end
	
	% If you threshold globally, what is the scale-dependent percentage
	% make up as a function of depth, out of the global number retained
	vwstats.tcount{index,jndex}=100*sum(~kilme(itshere))/nkept;      
      end % End loop over the scales
      
      % If not preconditioned, the individual scales are norm-preserving
      if prod(precon)==0
	difer(normscsq-vstats.normv^2,7,[],NaN)
	% Should perhaps report on how bad it gets - it's slight
	% But sure, this check gets done in the main body
      end
      
      % Also give the summary statistics for all of the levels
      vwstats.totnum{      J+1+1}=prod(size(vw));
      vwstats.percts{index,J+1+1}=prctile(vw(:),percs);
      vwstats.mjuvar{index,J+1+1}=[mean(vw(:)) var(vw(:))];
      vwstats.tcount{index,J+1+1}=100;
      % End of the extraction of scale-dependent information
      disp(sprintf('Finished layer %3.3i of %3.3i',index,ndeps))
    end
    % End of the loop over the layers, now save the results
    save(fname2,'vwstats','depkm','ndeps','percs','gnorjr','J','fs',...
	 'xlstring1','panel','wav','vw')
  end
  % End of the if that questions if calculations are needed
    
  % Check the size of the decomposition which shouldn't have budged
  difer(sum([vwstats.totnum{1:J+1}])-vwstats.totnum{J+2},[],[],NaN)
  % Check that the depths are the same; keep them in the struct though
  difer([vwstats.deptkm{:}]-depkm(:)',[],[],NaN)
  % Check that the proportions of the unthresholded coefficients adds up
  difer(sum(reshape([vwstats.tcount{:,1:J+1}],ndeps,[]),2)...
	-[vwstats.tcount{:,J+1+1}]',[],[],NaN)
  
  % What do the model norms add up to, which isn't one hundred?
  % This is NOT the partial correlation coefficient, which would remove
  % the effect of the other variables... 
  % http://faculty.vassar.edu/lowry/ch3a.html
  addup=sum(reshape([vwstats.normsq{:,1:J+1}],ndeps,[]),2);
  % This should be ABOUT hundred when ALL panels are considered as we are
  % comparing with the norm of ALL panels. If it's just ONE panel this is
  % the much more variable percentage of the contribution of that panel
  % to the overal structure
  if panel==0
    disp(sprintf('Norm squared percentages add up to %6.2f',...
		 mean(addup)))
  else
    disp(sprintf('Panel %i norm squared percentages add up to %6.2f',...
		 panel,mean(addup)))
  end

  % Only make the figure if there is no output
  if ~nargout
    % And now make the figure
    % The scale colors
    cols={'r','g','b','k','y',grey};
    % Which percentages will be PLOTTED
    percp=[85 100];
    [~,perci]=intersect(percs,percp);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf ; clear p
    ytix=[0 410 660 1000 2000 2899];
    ytil={'0' '410' '660' '1000' '2000' 'CMB'};
    
    [ah,ha]=krijetem(subnum(1,2));
    set(ah,'FontS',fs)
    axes(ah(1))
    % Shade the region which is the [variable] 85th to 100th percentile of
    % the entire set considered to provide a background for the percentiles
    % within the scales 
    allpercs=cat(1,vwstats.percts{:,J+1+1});
    lft=allpercs(:,1);
    rgt=allpercs(:,end);
    % Make sure it goes beyond the visible panel
    xtra=0.1;
    patchdep=[ytix(1)-xtra ; depkm ; ytix(end)+xtra];
    patchlft=[lft(1) ; lft ; lft(end)];
    patchrgt=[rgt(1) ; rgt ; rgt(end)];
    pp=patch([patchlft ; flipud(patchrgt)],...
	     [patchdep ; flipud(patchdep)],grey);
    delete(pp)
    hold on
    
    % Now plot the curves for the individual scales
    for jndex=1:J+1
      % The various percentiles
      allpercs=cat(1,vwstats.percts{:,jndex});
      for pndex=1:length(percp)
	p(pndex,jndex)=plot(allpercs(:,perci(pndex)),depkm,...
			    'Color',cols{jndex});
      end
      legs{jndex}=sprintf('wavs %i',jndex);
    end
    % xl=xlim; xlim([-1 xl(2)]); xlim([-1 25]); 
    xl=xlim; xlim([-1 xl(2)]); 
    if gnorjr==1
      xlim([-1 40]); 
    elseif gnorjr==2
      xlim([-1 80]); 
    end
    g{1}=plot(xlim,[ytix(2:end-1)' ytix(2:end-1)'],'k:');
    % Delete all but the 100th percentiles of the scales but keep all
    % percentiles of the total set
    set(p,'linew',1)
    hold off; delete(p(1:end-1,1:end))
    xl(1)=xlabel(xlstring1);
    yl(1)=ylabel('depth (km)');
    box on
    leg(1)=legend(p(2,:),legs,'Location','SouthEast');
    
    axes(ah(2))
    for jndex=1:J+1
      % The various percentiles
      alltcount=cat(1,vwstats.tcount{:,jndex});
      % Nope - the various rms proportions
      alltcount=cat(1,vwstats.normsq{:,jndex});
      p(pndex+1,jndex)=plot(alltcount,depkm,'Color',cols{jndex});
      hold on
      legs{jndex}=sprintf('scale %i',jndex);
    end
    % xl=xlim; xlim([0 xl(2)]); xlim([-2 50])
    xl=xlim; xlim([0 xl(2)]); xlim([-2 100])
    if ~~panel
      xlim([-2 30])
    end
    set(p(pndex+1,:),'linew',1)
    g{2}=plot(xlim,[ytix(2:end-1)' ytix(2:end-1)'],'k:');
    hold off
    % delete(p(end,end))
    % xlstring2=sprintf('%s contribution after %i%s thresholding',...
    % 		    '%',tperc,'%');
    xlstring2=sprintf('%s contribution to l2 model norm','\%');
    xlstring2='\textsf{contribution to $\ell_2$ model norm}';
    xl(2)=xlabel(xlstring2);
    set(xl(2),'Interpreter','LaTeX','FontS',fs+1)
    if ~~panel
      % Remember the reordering by the others
      pv=[5 4 6 2 3 1];
      legs=sprintf('chunk %i',pv(panel));
    end
    leg(2)=legend(legs,'Location','SouthEast');
    
    set(ah,'ydir','rev')
    set(ah,'ytick',ytix,'ytickl',ytil,'ylim',minmax(ytix))
    if ~~panel
      set(ah(2),'xtick',[0:10:30])
    else
      set(ah(2),'xtick',[0:25:100])
    end
    % Put on the letter labels
    [bh,th]=label(ha,'ur',10,2*[gnorjr==2]); 
    
    nolabels(ah(2),2)
    longticks(ah,2)
    fig2print(gcf,'portrait')
    shrink(ah,1.25,1)
    serre(ah,0.75,'across')
    movev(bh,2)
    movev(leg(1),.1)
    % Next line only for last modification
    movev(leg(1),.1)
    % Next line only for last modification
    movev(leg(2),.15)
    
    % Last minute
    if ~panel
      delete(leg(2))
    end
    set(findobj('string',sprintf('wavs %i',J+1)),...
	'string',sprintf('scals %i',J)) 

    % Identify panel
    if ~~panel
      wav=sprintf('%s_%i',wav,panel);
    end
    actprint=1;
    if gnorjr==1
      figna=figdisp([],sprintf('GN_%s',wav),[],actprint);
    elseif gnorjr==2
      figna=figdisp([],sprintf('JR_%s',wav),[],actprint);
    end
    if actprint==1
      system(sprintf('epstopdf %s.eps',figna));
    end
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Making sense? Check again
  if xver==1
    watdepth=ceil(randi(ndeps));
    watscale=J+1;
    allpercs=cat(1,vwstats.percts{:,watscale});
    al=log10(vw(:));
    hist(al(~isinf(al)),100)
    lk=log10(allpercs(watdepth,perci));
    xlim([-10 5])
    set(gca,'xtick',lk(~isinf(lk))); longticks
    grid on; rottick; shrink(gca,1.35,1.35)
  end
  % Potential output (it's vWstats!)
  % Let the first output be the panel number
  varns={panel,depkm,vwstats};
else
  % We're doing a single depth in this section and saving the wavelet
  % transforms for the next time. Don't forget that this is a recursive
  % portion of the algorithm. 
  if gnorjr==1
    % And now define what will come out of the results at the end
    fname=fullfile(getenv('IFILES'),'EARTHMODELS','MONTELLI','LORIS',...
		   sprintf('loris5_GN_%4.4i_%i_%i_%s_%i_%i.mat',...
			   round(depkm),N,J,wav,precon(1),precon(2)));
  elseif gnorjr==2
    % And now define what will come out of the results at the end
    fname=fullfile(getenv('IFILES'),'EARTHMODELS','RITSEMA','LORIS',...
		   sprintf('loris5_JR_%4.4i_%i_%i_%s_%i_%i.mat',...
			   round(depkm),N,J,wav,precon(1),precon(2)));
  end
  
  if exist(fname,'file')==2
    % Loading the wavelet transformed seismic models
    load(fname)
    if xver==1
      disp(sprintf('Loading %s',fname))
    end
  else
    if gnorjr==1
      disp('Reading Montelli model')
      % Load the model at a particular depth, faces rearranged and all
      v=readGNmodel(depkm);
    elseif gnorjr==2
      disp('Reading RITSEMA model')
      % Load the model at the requested, faces rearranged and all
      v=plm2cube(interpJRmodel(depkm));
      
    end
    % This should be on pixel centered, non-overlapping, completely covering,
    % registration with an even number of points. Right now we have this on
    % grid node registration, with an odd number of points. Later we will
    % change this, here we will fake this for the moment, to be quick
    disp('Temporary faking of the pixel center registration')
    v=v(1:2^N,1:2^N,:);
    
    % Save the norm etc, as when preconditioned, it changes
    % But remember this is taking ALL panels
    vstats.normv=norm(v(:));
    vstats.minmx=minmax(v(:));
    vstats.meanz=mean(v(:));
    vstats.variz=var(v(:));
    
    % Now calculate the wavelet transforms
    switch wav
     case 'D2'
      % The 2-tap Daubechies wavelet transform
      vw=angularD2WT(v,[J J],'forward',1);
     case 'D4'
      % The 4-tap Daubechies wavelet transform
      vw=angularD4WT(v,[J J],precon,'forward',1);
     case 'D6'
      % The 6-tap Daubechies wavelet transform
      vw=angularD6WT(v,[J J],precon,'forward',1);
    end
    
    % Quick plot?
    if xver==1
      clf
      [cb,xcb]=quickplot(N,J,vw,depkm);
    end
    
    % When you're done, save the results
    save(fname,'vw','depkm','wav','N','J','precon','vstats')
  end
  % Potential output (it's vstats!)
  varns={vw,depkm,vstats};
end

% Make output
varargout=varns(1:nargout);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cb,xcb,pgw]=quickplot(N,J,vw,depkm,dax)
% Font size of the label
defval('fs',8)
% Color map and saturation as percentiles
defval('colmap','kelicol');
% Grid information
wg.N=N; wg.J=J;
% Not at ALL the same as using HALVERANGE
colperc=[1 99];
% May go with HALVERANGE after all
defval('dax',round(halverange(vw,75,NaN)));
% Explicit and absolute color limits of the VALUE of the coeffs
defval('dax',prctile(vw(:),colperc));

% The actual plotting
[~,~,~,pgw]=plotoncube(vw,'2D',1,[],[],[],dax,[],0,100,wg);
% Color bar etc
colpos=[0.5616    0.1714+0.025    0.3143    0.0298];
[cb,xcb]=addcb(colpos,dax,dax,colmap,range(dax)/4);
set(cb,'fonts',fs)
set(xcb,'string',sprintf(...
    'wavelet coefficients at %i km',round(depkm)),'fonts',fs)

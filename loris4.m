function loris4(depkm,wav,N,J,precon,xver,doneven,actprint,gnorjr,agu)
% LORIS4(depkm,wav,N,J,precon,xver,doneven,actprint,gnorjr,agu)
%
% Map view (calculated on the fly) and curves illustrating
% quality-of-thresholded-reconstruction (calculated and saved and then
% loaded) showing the wavelet sparsity of an Earth model. 
%
% INPUT:
%
% depkm    The depth in km where the model is read; must match exactly
% wav      The chosen wavelet construction used, 'D2', 'D4' or 'D6'
% N        The power of the dyadic subdivision for the cubed sphere 
% J        Maximum scale in both directions for the wavelet analysis 
% precon   [1|0 1|0] toggles preconditioning in either dimension
% xver     1 Extra verification step and TOPOGRAPHY is being analyzed
%          0 Assuming that it all works, analyze MONTELLI's wave model
% doneven  0 It does all of the plotting and makes a nice figure
%          1 It only calculates the numbers and turns out a table
% actprint 1 Physically creates the plot and saves them to file
%          0 Suggests a suitable print command on the command line
% gnorjr   1 For the model by Montelli et al. (2006)
%          2 For the model by Ritsema et al. (2010)
% agu      1 For presentation-style quality
%          2 For SPIE
%
% EXAMPLES:
%
% loris4('demo1') % Makes something that can go in a booklet
% loris4('demo2') % Makes a table that can go in a paper
% loris4([],[],[],[],[],1) % Makes a graph of topography
% 
% actprint=0;
% loris4(406.350,'D4',7,3,[1 1],0,0,actprint,1) % Paper figure A
% loris4(406.350,'D4',7,3,[1 1],0,0,actprint,2) % Paper figure B
% loris4(677.250,'D4',7,3,[1 1],0,0,actprint,1,2) % SPIE figure 3A
% loris4(677.250,'D4',7,3,[1 1],0,0,actprint,2,2) % SPIE figure 3B
%
% Last modified by fjsimons-at-alum.mit.edu, 8/29/2011

% Identify defaults
defval('depkm',406.350);
defval('depkm',1015.875);
defval('depkm',609.525);
% If 1, GN model, if 2, JR model
defval('gnorjr',1)
defval('agu',0)

if ~isstr(depkm)
  defval('wav','D4')
  defval('N',7)
  defval('J',N-4)
  defval('precon',[1 1]);
  defval('xver',0)
  defval('doneven',0)
  defval('actprint',0);

  % This for the legend
  prec={'raw','precond'};

  % The thresholding percentages that we will be plotting as MAPS
  tpercs1=[0 50 85 95];
  if agu==2
    % For SPIE, make this fewer
    tpercs1=[0 85 95];
  end
  
  % How many percentiles will you ultimately calculate?
  lst=100;
  
  if xver==0
    if gnorjr==1
      disp('Reading Montelli model')
      % Load the model at the requested, faces rearranged and all
      v=readGNmodel(depkm);
      % This is where the output goeth
      dirname=fullfile(getenv('IFILES'),'EARTHMODELS','MONTELLI','LORIS'); 
      % And now define what will come out of the results at the end
      fname=fullfile(dirname,sprintf('loris4_GN_%4.4i_%i_%i_%s_%i_%i_%i.mat',...
				     round(depkm),N,J,wav,precon(1),precon(2),lst));
      cbstring=sprintf('%s wave speed anomaly (%s)','P','%');
      upv=1.5;
      dnv=-1.5;
      if agu==2
	% In SPIE
	upv=1.25;
	dnv=-1.25;
      end
    elseif gnorjr==2
      disp('Reading RITSEMA model')
      % Load the model at the requested, faces rearranged and all
      v=plm2cube(interpJRmodel(depkm));
      % This is where the output goeth
      dirname=fullfile(getenv('IFILES'),'EARTHMODELS','RITSEMA','LORIS');
      % And now define what will come out of the results at the end
      fname=fullfile(dirname,sprintf('loris4_JR_%4.4i_%i_%i_%s_%i_%i_%i.mat',...
				     round(depkm),N,J,wav,precon(1),precon(2),lst));
      cbstring=sprintf('%s wave speed anomaly (%s)','S','%');
      upv=2.5;
      dnv=-2.5;
      if agu==2
	% In SPIE
	upv=2;
	dnv=-2;
      end
    end
  else
    % In this extra-verification step we work with topography to make sure
    % we are seeing and doing things right
    defval('L',ceil(2^(N+1)))
    % Yes, that should be loris2 in this filename, see LORIS2
    fname=fullfile(dirname,sprintf('loris2_%i_%i.mat',N,L));
    load(fname)
    % And now define what will come out of the results at the end
    fname=fullfile(dirname,sprintf('loris4_%i_%i_%i_%s_%i_%i_%i.mat',...
			   N,L,J,wav,precon(1),precon(2),lst));
  end

  % This should be on pixel centered, non-overlapping, completely covering,
  % registration with an even number of points. Right now we have this on
  % grid node registration, with an odd number of points. Later we will
  % change this, here we will fake this for the moment, to be quick
  if xver==1 
    disp('Temporary faking of the pixel center registration');
  end
  v=v(1:2^N,1:2^N,:);

  if doneven==0
    % Start the figure
    clf
    [ah,ha]=krijetem(subnum(length(tpercs1)+1,2));
  end

  % Do the calculations that result in MAPS being plotted
  for index=1:length(tpercs1)
    % Perform the thresholded reconstruction on all faces
    [bl,vw,vwt,recer1(index,:),kilme]=...
	angularthresh(v,tpercs1(index),J,0,wav,precon,0);

    % Make messages for what you're plotting or at the least displaying
    mes1=sprintf('%2i%s thresholded\n %2i%s error norm',...
		 tpercs1(index),'%',round(recer1(index,1)),'%');
    mesa=sprintf('%2i%s & %6.3f',...
		 tpercs1(index),'%',recer1(index,1));
    mes2=sprintf('%s %s\n %i nonzero',...
		 wav,prec{prod(precon)+1},sum(~kilme(:)));
    mesb=sprintf('%s %s %i nonzero',...
		 wav,prec{prod(precon)+1},sum(~kilme(:)));
    mes3=sprintf('%4i km\n N = %i J = %i',round(depkm),N,J);
    mesc=sprintf('%4i km N = %i J = %i',round(depkm),N,J);
    
    if doneven==0
      % Project to Mollweide - this changes the values!
      [lonm,latm,vd]=cube2moll(bl);
      
      % Adjust the color range - after the projection!
      if xver==0
	disp(sprintf('Old model range %4.1f to %4.1f',min(vd(:)),max(vd(:))))
	vd(vd>upv)=upv;
	vd(vd<dnv)=dnv;
	disp(sprintf('New model range %4.1f to %4.1f',min(vd(:)),max(vd(:))))
      end
      
      % Make a map
      axes(ha(index))
      pc=pcolor(lonm,latm,vd); shading flat; 
      
      % Put a color bar for the tomography 
      if xver==0 && index==1
	pos=[0.2042 0.7765 0.1875 0.0108];
	cb=colorbarf('hor',10,'Helvetica',pos);
	set(cb,'xtick',linspace(dnv,upv,5),'xlim',[dnv upv])
	longticks(cb)
	axes(cb)
	xlcb=xlabel(cbstring);
	axes(ha(index))
        % Note get(gcf,'NextPlot')
      end

      if xver==0
	kelicol
      else
	% Plot in the two-dimensional plane - no reordering needed
	[dem,dax,ziro]=sergeicol;
	% Hardwire so that the zero crossing is visually pleasing
	dax=[-7473 5731].*[1.025 0.975];
	caxis(dax)
	colormap(sergeicol)
      end
      
      % These are the continents and plate boundaries in mollweide
      [a,b]=plotcont([],[],2); 
      [c,d]=plotplates([],[],2); 
      % plotoncube(v,'2D',[],[],[],[],[],[],0); plotcont([],[],9)
      axis image off
      axis(axis*1.05)
      tx (index)=text( 2.25,-1.25,mes1);
      tx2(index)=text(-2.25,-1.25,mes2);
      if xver==0
	tx3(index)=text(2.25,1.25,mes3);
      end
    else
      if index==1; 
	disp(mesb); disp(mesc); disp(' '); 
      end
      disp(mesa)
    end
  end
  
  if doneven==0
    % Now move on to plotting the curve
    if exist(fname,'file')==2
      load(fname)
      disp(sprintf('%s loaded',fname))
    else
      % The thresholding percentages that we will be calculating
      tpercs2=linspace(0,99,lst);
      
      % Start the Matlab CPU pool if you have more than one processor
      % Then change the for to parfor. However, are there savings?
      matlabpool open
      
      % Start the diagnostics vector
      recer2=nan(length(tpercs2),6);
      
      % Do the calculations that result in the GRAPH being plotted
      tic
      for index=1:length(tpercs2) 
	disp(sprintf('Performing wavelet transforms %3.3i / %i',...
		     index,length(tpercs2)))
	% Perform the thresholded reconstruction on all faces
	[~,~,~,recer2(index,:)]=...
	    angularthresh(v,tpercs2(index),J,0,wav,precon,0);
      end
      toc
      
      % Make sure that if we are NOT requesting truncation, we get perfect
      % reconstruction and all other stuff makes sense
      difer(recer2(find(~tpercs2),[1 2 4]),9)
      
      matlabpool close
      if xver==0
	% Save the results for future use
	save(fname,'recer2','tpercs2','N','J','wav','lst','precon')
      else
	% Save the results for future use
	save(fname,'recer2','tpercs2','N','L','J','wav','lst','precon')
      end    
      disp(sprintf('%s saved',fname))
    end
    
    % Now plot the diagnostics in the final panel
    final1=length(tpercs1)+1;
    axes(ha(final1))
    
    % Instead of the percentile of the data which we used to truncate, now
    % use the actual percentage of the coefficients that is being
    % zeroed. Unless there are repeated values, which is unlikely, these two
    % numbers are expected to be identical, but nevertheless 
    tpercs3=recer2(:,2);
    
    % Actually, the plot will look better without the 0 [and 100] marks on it
    % is it is to make it to log scale
    % recer2=recer2(find(tpercs3),:);
    % tpercs3=tpercs3(find(tpercs3),:);
    
    % Plot the error of the reconstruction, in percent
    p(1)=plot(tpercs3,recer2(:,1),'k-');
    % The offset controls how far down the last label box needs to go
    ofs=[0 0 0 1.5+2*[max(recer2(:,6))>101]];
    lims1=[-3 102 -1 25];
    axis(lims1-ofs)
    set(ha(final1),'ytick',[0:5:25])
    hold on
    % Add eplicit grid lines
    plot(xlim,[0 0],':')
    plot([tpercs1(:) tpercs1(:)],ylim,':')
    p(2)=plot(tpercs1,recer1(:,1),'ko');
    hold off
    xl(1)=xlabel('percentile thresholded (%)');
    yl(1)=ylabel('l2 error norm (%)');
    shrink([ha(final1)],1.1,1)
    % Now add the second axis
    axi(1)=xtraxis(ha(final1),tpercs1,tpercs1);
    % Plot the L1-NORM of the wavelet coefficients
    p(3)=plot(tpercs3,recer2(:,5),'b');
    hold on
    p(4)=plot(tpercs1,recer1(:,5),'bo');
    % Plot the TOTAL VARIATION on the same axis
    p(5)=plot(tpercs3,recer2(:,6),'r');
    p(6)=plot(tpercs1,recer1(:,6),'ro');
    hold off
    set(axi(1),'yaxisl','r')
    lims2=[lims1(1:2) 50 101.5+5*[max(recer2(:,6))>101]];
    axis(lims2)
    set(axi(1),'Color','none')
    nolabels(axi(1),1)
    noticks(axi(1),1)
    yl(2)=ylabel('l1 norm | total variation (%)');
    set(yl(2),'rotation',-90)
    moveh(yl(2),5)
    set(axi(1),'box','off','xaxisl','t')
    
    % Cosmetics etc
    if xver==0
      serre(ha(2:final1-1),1/3,'down')
    end
    longticks([ha(final1) axi(1)],1)
    set(tx,'horizontala','left')
    set(tx2,'horizontala','right')
    set(ha(final1),'xtick',[0:20:80 99])
    set(p(2),'MarkerF','k','MarkerS',4)
    set(p(4),'MarkerF','b','MarkerS',4)
    set(p(6),'MarkerF','r','MarkerS',4)
    fig2print(gcf,'tall')
    % Put on the letter labels
    [bh,th]=label(ha,'ul',6,[length(tpercs1)+1]*[gnorjr==2]); 
    set(th,'FontS',20)
    % But now lower the label effectively
    set(ha(final1),'ylim',lims1(3:4))
    delete(ha(final1+1:end))
    axes(axi(1))
    
    % Do the printing to EPS
    popts='-zbuffer'',''-r300';
    if xver==1
      figna=figdisp([],wav,popts,actprint);
    else
      if gnorjr==1
	figna=figdisp([],...
		      sprintf('GN_%s_%4.4i',wav,round(depkm)),popts, ...
		      actprint);
      elseif gnorjr==2
	figna=figdisp([],...
		      sprintf('JR_%s_%4.4i',wav,round(depkm)),popts, ...
		      actprint);
      end
    end
    if actprint==1
      % Convert to PDF
      system(sprintf('epstopdf %s.eps',figna));
    end
  end

  if agu==2
    movev(cb,-.025)
    movev(ah(3),.02)
    movev(ah(5),.04)
    movev(ah(7),.08)
    movev(axi(1),.08)
    shrink([ah(7) axi(1)],1,1.25)
    movev([cb ah([1 3 5 7]) axi(1)],.05)
  end

  if agu==1
    fig2print(gcf,'landscape')
    movev(ha(2),0.205)
    moveh(ha(2),0.4)
    movev(ha(3),0.15)
    movev(ha(4),0.15*2+0.005)
    moveh(ha(4),0.4)
    movev(ha([3 4]),-.02)
    movev([ha(5) axi(1)],0.15*2-0.05)
    moveh([ha(5) axi(1)],0.2)
    set(th(1:5),'FontS',12)
    delete([bh(1:5) th(1:5)])
    % Reduce more folks
    % delete(ha(3:4))
    delete(ha(2:3))
    movev(ha(4),0.205+0.025+0.002)
    movev([ha(5) axi(1)],0.15)
    % shrink([ha(5) axi(1)],1.5,1)
    %set([ha(1:2)],'CameraV',3)
    set([ha([1 4])],'CameraV',3)
    movev(cb,-.05)
    shrink([ha(5) axi(1)],1,0.75)
    moveh([ha(1) cb],-.0125)
    %moveh(ha(2),.025)
    moveh(ha(4),.025)
    
    if gnorjr==1
      figna=figdisp([],...
		    sprintf('GN_%s_%4.4i_agu',wav,round(depkm)),popts, ...
		    actprint);
    elseif gnorjr==2
      figna=figdisp([],...
		    sprintf('JR_%s_%4.4i_agu',wav,round(depkm)),popts, ...
		    actprint);
    end
    if actprint==1
      system(sprintf('epstopdf %s.eps',figna));
    end
  end
elseif strcmp(depkm,'demo1')
  % d=readGNmodel(NaN);  
  d=load(fullfile(getenv('IFILES'),'EARTHMODELS','MONTELLI','GNdepths'));
  % In km now
  d=d/1000;
  % Ready for merging into a pdf booklet via mergepdf
  for i=[1 3 6 10 15 19 28 31 46 68 90 129];
    loris4(d(i),'D2',[],[],[],[],[],1); 
    loris4(d(i),'D4',[],[],[],[],[],1); 
    loris4(d(i),'D6',[],[],[],[],[],1); 
  end
elseif strcmp(depkm,'demo2')
  % d=readGNmodel(NaN);  
  d=load(fullfile(getenv('IFILES'),'EARTHMODELS','MONTELLI','GNdepths'));
  % In km now
  d=d/1000;
  % Ready for tabulating some special values
  % Ritsema model
  for i=[10 19 28 46 90];
    disp(' ')
    loris4(d(i),'D2',[],[],[],[],1,[],2); 
    disp(' ')
    loris4(d(i),'D4',[],[],[],[],1,[],2); 
    disp(' ')
    loris4(d(i),'D6',[],[],[],[],1,[],2); 
    disp('---------------------------------------')
  end
  % Ready for tabulating some special values
  % Montelli model
  for i=[10 19 28 46 90];
    disp(' ')
    loris4(d(i),'D2',[],[],[],[],1,[],1); 
    disp(' ')
    loris4(d(i),'D4',[],[],[],[],1,[],1); 
    disp(' ')
    loris4(d(i),'D6',[],[],[],[],1,[],1); 
    disp('---------------------------------------')
  end
end

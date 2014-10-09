function treerings(diro,fname)
% TREERINGS(diro,fame)
%
% Analizes a treering scan. 
%
% INPUT:
%
% diro       A directory name string, e.g. '.'
% fname      A filename string, e.g. 'KK3_barkup.tif'
%
% EXAMPLE:
%
% treerings('.',[]) % Current directory - default file
% treerings([],' KJ10_barkup.tif') % Default directory, some file
% treerings('/Volumes/home/fjsimons/','somefile') % You get the idea 
%
% Last modified by fjsimons-at-alum.mit.edu, 12/12/2008

% Where are they kept? This will be different for you
defval('diro','/home/fjsimons/CLASSES/FRS145/2008/FieldData/TreeRings');

% Which filename is it? I supply a default - you can delete this line
defval('fname','SJ2_barkup.tif')
% Do the analysis or not
defval('anal',1)
% Number of x tick marks
defval('nticks',10)
% Number of sections to display for analysis
defval('nsex',8)
% Markersize in section plots
defval('marks',6)
% Markersize in overview plot
defval('marksa',6)
% Actually make all the plots, or just pretend
defval('plotit',1)

% Read in the tree scan...
rgb=imread(fullfile(diro,fname));

% Make sure the long side is across the screen
if size(rgb,1)>size(rgb,2)
  % Transpose...
  red=rgb(:,:,1)';
  green=rgb(:,:,2)';
  blue=rgb(:,:,3)';
  rgb=zeros([size(red) 3],'uint8');
  rgb(:,:,1)=red;
  rgb(:,:,2)=green;
  rgb(:,:,3)=blue;
else
  % Or don't...
  red=rgb(:,:,1);
  green=rgb(:,:,2);
  blue=rgb(:,:,3);
end

% These values are unsigned 8-bit integer (from 0-255) so they require
% somewhat special attention to convert to grey scales
grae=uint8([double(red)+double(green)+double(blue)]/3);

% Select the MIDPOINT on the y-axis for plotting and guides
ypoint=round(size(grae,1)/2);

clf
% Plot the color image
ah(1)=subplot(4,1,1);
ob(1)=image(rgb);
yl(1)=ylabel(sprintf('color'));
tl(1)=title(nounder(fname));
set(tl,'FontS',15)
% Plot the guide line
hold on
pg(1)=plot(xlim,[ypoint ypoint],'k');
hold off

% Plot the gray-scale image
ah(2)=subplot(4,1,2);
ob(2)=image(repmat(grae,[1 1 3]));
yl(2)=ylabel(sprintf('grey'));
% Plot the guide line
hold on
pg(1)=plot(xlim,[ypoint ypoint],'k');
hold off

% Plot the derivative of the gray-scales as a time series
ah(3)=subplot(4,1,3);
deriv=diff(double(grae(ypoint,:)));
ob(3)=plot(abs(deriv),'k');
axis tight
% Saturate this axis to show the largest changes - you may change this
ylim([prctile(abs(deriv),75) max(abs(deriv))])
yl(3)=ylabel(sprintf('%sgrey @ %i','\Delta',ypoint));

% Plot the gray-scales as a time series
ah(4)=subplot(4,1,4);
ob(4)=plot(grae(ypoint,:),'k');
axis tight
yl(4)=ylabel(sprintf('grey @ %i',ypoint));
xels=xlim;
set(ah(4),'xtick',round(linspace(xels(1),xels(2),nticks)))
xl=xlabel('pixel number');
ymean=mean(grae(ypoint,:));
hold on
plot(xels,[ymean ymean],':')
hold off

% Some cosmetic changes
set(ah,'ytick',[])
nolabels(ah(1:3),1)
set(ah(3:4),'xgrid','on','ygrid','on')
longticks(ah,4)
fig2print(gcf,'landscape')
% Print figure
figdisp(nounder(pref(fname)),[],[],plotit)

% Now do the analysis
if anal==1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  r=input(sprintf('\nReady for analysis? [1/0]\n'));
  if r==1
    disp(sprintf('\nClick on the figure ; end series with return\n'));
    disp(sprintf('x pixel location and mouse button are recorded\n'))
    
    % Look at the figure quarter by quarter
    xells=pauli(linspace(xels(1),xels(2),nsex+1))+...
	 [0 0 ; ones(nsex-1,1) zeros(nsex-1,1)];
    % And make the picks
    for index=1:nsex
      set(ah,'xlim',xells(index,:),'xtick',...
	     round(linspace(xells(index,1),xells(index,2),nticks)))
      % Adjust title
      axes(ah(1))
      tl(1)=title(sprintf('%s section %i / %i',...
			  nounder(fname),index,nsex));
      % Plot new mean of the grey values
      axes(ah(4))
      hold on
      plot(xells(index,:),mean(grae(ypoint,...
			    round(xells(index,1)):round(xells(index,2)))),':')
      hold off
      
      % Saves x and y coordinate, and mouse button pressed
      [xpicks{index},ypicks{index},bpicks{index}]=ginput;
    end
    % And then collect all the picks together - to the nearest pixel
    allpicks=round(cat(1,xpicks{:}));
    ballpicks=round(cat(1,bpicks{:}));
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  r=input(sprintf('\nReady for section review? [1/0]\n'));
  if r==1
    % Review what you've picked
    % Different symbols for the mousebuttons
    symbs={'v','s','^'};
    for index=1:nsex
      set(ah,'xlim',xells(index,:),'xtick',...
	     round(linspace(xells(index,1),xells(index,2),nticks)))
      axes(ah(1)); hold on
      tl(1)=title(sprintf('%s section %i / %i',...
			  nounder(fname),index,nsex));
      pp{1}=repmat(NaN,1,length(xpicks{index}));
      for ondex=1:length(xpicks{index})
	pp{1}(ondex)=plot(xpicks{index}(ondex),...
			  ypoint,symbs{bpicks{index}(ondex)});
      end
      hold off
      
      axes(ah(2)); hold on
      pp{2}=repmat(NaN,1,length(xpicks{index}));
      for ondex=1:length(xpicks{index})
	pp{2}(ondex)=plot(xpicks{index}(ondex),...
			  ypoint,symbs{bpicks{index}(ondex)});
      end
      hold off
      
      % Cosmetics
      set(pp{1},'MarkerF','w','MarkerE','k','MarkerS',marks)
      set(pp{2},'MarkerF','w','MarkerE','k','MarkerS',marks)
      
      axes(ah(3))
      set(ah(3),'xtick',unique(round(xpicks{index})))
      
      axes(ah(4))
      nolabels(ah(4),1)
      set(ah(4),'xtick',unique(round(xpicks{index})))

      % Any more picks? While you're at it, keep adding
      r=input(sprintf('\nAny more picks to be added? [1/0]\n'));
      if r==1
	[xadd{index},yadd{index},badd{index}]=ginput;
	axes(ah(1)); hold on
	ppadd{1}=repmat(NaN,1,length(xadd{index}));
	for ondex=1:length(xadd{index})
	  ppadd{1}(ondex)=plot(xadd{index}(ondex),...
			    ypoint,symbs{badd{index}(ondex)});
	end
	hold off
	axes(ah(2)); hold on
	ppadd{2}=repmat(NaN,1,length(xadd{index}));
	for ondex=1:length(xadd{index})
	  ppadd{2}(ondex)=plot(xadd{index}(ondex),...
			    ypoint,symbs{badd{index}(ondex)});
	end
	hold off

	% Cosmetics
	set(ppadd{1},'MarkerF','w','MarkerE','k','MarkerS',marks)
	set(ppadd{2},'MarkerF','w','MarkerE','k','MarkerS',marks)
	
	axes(ah(3))
	set(ah(3),'xtick',unique(round([xpicks{index} ; xadd{index}])))
	
	axes(ah(4))
	nolabels(ah(4),1)
	set(ah(4),'xtick',unique(round([xpicks{index} ; xadd{index}])))
      
	% And then collect all the picks together - to the nearest pixel
	allpicks=[allpicks ; round(xadd{index})];
	ballpicks=[ballpicks ; round(badd{index})];
	pause
      end
      
      % Print section figures with picks
      figdisp(nounder(pref(fname)),sprintf('%i_pix',index),[],plotit)
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  r=input(sprintf('\nReady for final review? [1/0]\n'));
  if r==1
    % Change size of markers
    set(findobj('MarkerS',marks),'MarkerS',marksa)
    set(ah,'xlim',xels,'xtick',...
	   round(linspace(xels(1),xels(2),nticks)))
    axes(ah(1))
    tl(1)=title(nounder(fname));
    % Print final figure with all picks
    figdisp(nounder(pref(fname)),'pix',[],plotit)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %r=input(sprintf('\nReady for data save? [1/0]\n'));
  %if r==1
    % APPEND these button and x-only data picks to an ascii file
    fid=fopen(sprintf('%s.pix',nounder(pref(fname))),'w');
    % Sort all of the picks by x - don't forget to also sort the buttons 
    [allpicks,bi]=sort(allpicks);
    a=fprintf(fid,'%i %i\n',[ballpicks(bi) allpicks]');
    fclose(fid);
  %end
end

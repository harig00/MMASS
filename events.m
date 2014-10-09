function events(ddir,dec,wsec,wolap,twoplot)
% EVENTS(ddir,dec,wsec,wolap,twoplot)
%
% Makes diagnostic plots of pairs of seismograms contained in a
% directory. This is a front end that passes filenames and options to
% the more useful function SIGNALS.
%
% INPUT:
%
% ddir      Directory in which the data. If the directory name contains: 
%           EVENTS or HRSECTIONS, special actions are taken.
% dec       0 Nothing special (default)
%           1 Reads decimated data from file *_dec.sac
%           2 Data in file will be low-passed on the fly
% wsec      Window length in second [default: 5]
% wolap     Fraction of window overlap [default: 0.875]
% twoplot   1 If you have the entire data set, will also plot this 
%           0 Only one signal plot is being generated [default]
%
% SEE ALSO:
%
% SIGNALS, DATANAL1, DATANAL2, DATANAL3
%
% Last modified by fjsimons-at-alum.mit.edu, 05/24/2010

defval('ddir','/home/fjsimons/MERMAID/09102004/EVENTS/')
defval('ddir','/home/fjsimons/MERMAID/11042003/EVENTS/')
defval('ddir','/home/fjsimons/MERMAID/08092007/EVENTS/')
defval('dec',0)
defval('wsec',5)
defval('wolap',0.875);
defval('twoplot',0);

% Divide the entire directory listing into file pairs 
names=filepair(ddir);

if twoplot==1
  figure(1); clf
  if strfind(ddir,'11032003')
    % Plot entire first data series
    datanal1([],[],[],[],1,0,0);
    xlim([250 150000])
    ylim([-6 -3.5]*1e4)
  elseif strfind(ddir,'09102004')
    % Plot entire second data series
    datanal2([],[],[],[],1,0,0);
    xlim([400 146000])
    ylim([-9.5 -3]*1e3)
  elseif strfind(ddir,'08092007')
    % Plot entire third data series
    datanal3([],[],[],[],1,0);
    xlim([0 167400])
    ylim([-2 2.5]*1e5)
  end
  % Get rid of the ridiculous axis markings
  delete(findobj('color',[.7 .7 .7]))
end

% Attempt at cross-referencing with the catalog
if ~isempty(findstr(ddir,'EVENTS'))
  load(fullfile(ddir,'../CATALOGS/NEIC9'))
end

% Uses SIGNALS to do the plotting and everything else
for index=1:length(names)
  % Get the time of the first arrival from the filename
  tim1=findtime(names{index}(1),ddir);
  % Get the time of the second arrival from the filename
  if length(names{index})>1
    tim2=findtime(names{index}(2),ddir);
  end
  % Now plot these on the first plot if it exists
  if twoplot==1
    % Adaptation for two plots at a time
    figure(1)
    if exist('xp','var'); delete(xp); end
    xp(1)=plot([tim1 tim1],ylim,'k-');
    xp(2)=plot([tim2 tim2],ylim,'k--');
    figure(2); clf
  end

  % Now make the actual signal plot
  signals(0,dec,ddir,names{index},wsec,wolap,1,0);

  % Now we have the timing information for P, cross-reference this
  % with the distance and the magnitude and also plot this on there
  if exist('NEIC9')
    data1=NEIC9(NEIC9==tim1,:);
    data2=NEIC9(NEIC9==tim2,:);
    str=sprintf('M1 %4.2f %s %3.1f / M2 %4.2f %s %3.1f',...
		data1(2),'\Delta',data1(3),...
		data2(2),'\Delta',data2(3));
    figc
    a=text(1-0.9412,0.02,str,...
	   'vert','bottom','horiz','left','fontsize',9);
    axisc
    figdisp(sprintf('signals_%s_%i',pref(names{index}{1}),dec))
  end
  pause
end

% Auxiliary function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ztime=findtime(zname,zdir)

if isempty(cell2mat(strfind(zname,'_')))
  % Filename mpilot??????.sac
  ztime=str2num(cell2mat(pref(pref(zname),'mpilot')));
else
  % For the new filename convention
  % Filename mpilot?_??????.sac
  ztime=str2num(suf(pref(zname),'_'));
end
if ~isempty(findstr(zdir,'HRSECTIONS'))
  % One hour with 10 minutes overlap!
  ztime=ztime*3600-(index-1)*600; 
end

function wtALL
% WTALL
%
% Analyzes all events in a directory by a wavelet transform.
% Simply looks at the first significant wavelet coefficient in the
% interval supplied in the SAC files.
%
% Last modified by fjsimons-at-alum.mit.edu, 11/30/2006

% Identify the location of the SAC files
dirs='/home/fjsimons/EALARMS';
% This directory lists earthquake codes followed by one called 'BHZdata'
% in which the SAC files reside
cd(dirs)

% Get rid of initial two nonsensival directory listings
diro=indeks(ls2cell(fullfile(dirs,'*.*')),'3:end');

% If plots is 1, plots, no data written
plots=0;
% Thresholding level
tlev=1;

if plots==0
  fid=fopen('wtALL_diagnostics','a');
end

rind=0;
for index=1:length(diro)
  cd(fullfile(diro{index},'BHZdata'))
  disp(' ')
  disp(num2str(index))
  disp(' ')
  disp(fullfile(diro{index},'BHZdata'))
  seis=ls2cell('*.sac.*iso');
  for ondex=1:length(seis)
    % Here perhaps could investigate the difference between GCARC, DIST,
    % and my own estimate... without ellipticity 
    disp(seis{ondex})
    [s,h,jk1,jk2,t]=readsac(seis{ondex},0);
    % Record length. in seconds
    rlen=(h.NPTS-1)*h.DELTA;
    % FILTERING THE DATA TO WITHIN 10 HZ DOES NOT DO MUCH AT ALL
    % LOW-PASS NOR BANDPASS
    % DEMEAN THE DATA SO THE EDGE EFFECTS ARE MINIMIZED
    % CHECK THE PLOT IS RIGHT, COMPARE WITH WAVELETS1, WAVELETS2
    % This performs the wavelet transform using the [2,4] construction
    [a,d,an,dn]=wt(detrend(s,'constant'),'CDF',[2 4],5,4);
    % MAKE SURE THE LAST LEVEL HAS PLENTY OF COEFFICIENTS IN IT
    % AND THAT THAT "FOUND" ONE IS NOT A BOUNDARY COEFFICIENT
    % IN OTHER WORDS, NEED TO WATCH FOR EDGE EFFECTS HERE.
    
    % What is the support? Zero all but one at the middle and see for
    % yourself. 
    
    % THRESHOLD THE WAVELET COEFFICIENTS
    [dt,dnz]=threshold(d,dn,'soft',tlev);
    % Now find the relative coefficients... BLIND!
    % With complete disregard for the signal to noise ratio
    % First nonzero coefficients
    for undex=1:length(dt)
      foundit=min(find(abs(dt{undex})));
      if ~isempty(foundit)
	% Returns the ABSOLUTE magnitude of the wavelet coefficient!
	cof(undex)=abs(dt{undex}(foundit));
      else
	cof(undex)=NaN;
      end
    end
    % Cumulative seismogram number
    rind=rind+1;
    % Make a plot if you want... to verify
    if plots==1
      ah(1)=subplot(211);
      p(1)=plot(t,detrend(s,'constant'),'-'); axis tight
      hold on
      yli=ylim;
      pp(1)=plot([h.T0 h.T0],[yli(1) yli(2)],'k--');
      hold off
      ah(2)=subplot(212);
      [pdy,stp,rgp]=dyadplot(detrend(s,'constant'),a,dt,an,dn,1,[h.B h.E]);
      hold on
      pp(1)=plot([h.T0 h.T0],[yli(1) yli(2)],'k--');
      hold off
      colormap(flipud(gray(128)))
      pause
    else
      % Write the diagnostics to a file
      % How close is SAC's DIST to GCARC*fralmanac('DegDis')...?
      fprintf(fid,['%5i %4.2f %3i %5.2f %8.2f %5i'...
		   repmat('%12.3f',1,length(cof)) '\n'],...
	      rind,rlen,tlev,h.MAG,h.DIST,index,cof);
    end
  end
  cd('..')
  cd('..')
end
fclose(fid);

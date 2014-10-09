function varargout=grace2plmt(Pcenter,units,forcenew)
% [potcoffs,cal_errors,thedates,datanames,errornames]=GRACE2PLMT(Pcenter,units,forcenew)
%
% This program reads in the Level-2 GRACE geoid products from either the CSR or
% GFZ data centers, does some processing, and saves them as a plmt matrix
% in a .mat file.  In particular, the coefficients are reordered to our
% prefered lmcosi format, they are referenced to the WGS84 ellipsoid, and
% the C2,0 coefficients are replaced with more accurate measurements from
% satellite laser ranging.  You have the option of leaving them as geopotential 
% or converting them to surface mass density using the method of 
% Wahr et al. 1998, based on Love numbers (see PLM2POT).
%
% INPUT:
% 
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
% units       'POT' or 'SD' for whether you want geopotential or surface
%               mass density
% forcenew    Wether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
%
% OUTPUT:
% 
% Returns these variables and saves them in a .mat file:
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 6]
%                    these could also be in surface mass density
%    cal_errors     calibrated errors [nmonths x addmup(Ldata) x 4]
%    thedates       time stamps in Matlab time
%    datanames      filenames from which the potential coefficients came
%    errornames     filenames from which the calibrated errors came
%
% NOTE:
%   SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/J2/ e.g.
%   ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL04_2010_12.txt
%   ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL04_2012_03.txt
%  The header was removed and the file renamed to SLR_C20.txt for easy use.
%  Updated files keep getting posted in the same location.  
%
% Last modified by charig-at-princeton.edu, 06/26/2012
% Last modified by fjsimons-at-alum.mit.edu, 01/14/2013

% Determine parameters and set defaults
defval('Pcenter','CSR')
defval('units','POT')
defval('forcenew',0)

% Top level directory
% For Chris
IFILES=getenv('IFILES');
% For FJS, who has a different $IFILES
%IFILES='/u/charig/Data/';

% Where the original data files are kept
%defval('ddir1',fullfile(IFILES,'GRACE','Originals',Pcenter));
defval('ddir1',fullfile(IFILES,'GRACE',Pcenter));

% Where you would like to save the new .mat file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE'));
% And the name of that save file
if strcmp(units,'SD')
    fnpl=sprintf('%s/%s_alldata_%s.mat',ddir2,Pcenter,units);
else
    fnpl=sprintf('%s/%s_alldata.mat',ddir2,Pcenter);
end

% If this file already exists, load it.  Otherwise, or if we force it, make
% a new one (e.g. you added extra months to the database).
if exist(fnpl,'file')==2 && forcenew==0
     load(fnpl)
     disp(sprintf('%s loaded by GRACE2PLMT',fnpl))
else

% DATA CENTER
if Pcenter == 'GFZ'
   % Find the coefficient files
   datanames=ls2cell(fullfile(ddir1,'GSM*G---_0004'));
   % Find the error files
   errornames=ls2cell(fullfile(ddir1,'GSM*G---_0004.txt'));
   % Know a priori what the bandwidth of the coefficients is
   Ldata=120; 
elseif  Pcenter == 'CSR'
   datanames=ls2cell(fullfile(ddir1,'GSM*0060_0004'));
   errornames=ls2cell(fullfile(ddir1,'GSM*0060_0004.txt'));
   % Know a priori what the bandwidth of the coefficients is
   Ldata=60;
end

% WGS84 reference SETUP
% For now just hardcode the even zonal coefficients (J), later use
% Frederik's GRS.m program, don't bother with the higher degrees
j2= 0.108262982131e-2*-1.0/(2*2+1)^0.5; % will be row 4
j4=-0.237091120053e-5*-1.0/(2*4+1)^0.5; % will be row 11

% C20 CORRECTION SETUP

% Load the C(2,0) coefficients from satellite laser ranging
slrc20=load(fullfile(IFILES,'GRACE','SLR_C20.txt'));
% The sigma error is column 4
slrc20_error=slrc20(:,4)*1e-10;
% Remove the AOD1B model which was removed from the GRACE GSM data but
% restored to the SLR data.  Use the raw value (column 2). See webpage.
slrc20=[slrc20(:,1) slrc20(:,2)-slrc20(:,5)*1e-10];
% Convert the dates to Matlab format
[n,m]=size(slrc20);
slrc20(:,1)=datenum([slrc20(:,1) ones(n,1) ones(n,1)]);
% Make slrc20 relative to the WGS84 ellipsoid
slrc20(:,2) = slrc20(:,2) - j2;

% Initialize
nmonths = length(datanames);
thedates = zeros(1,nmonths);
[dems,dels]=addmon(Ldata);

% Calibrated errors are normally used instead, but they are kept here for
% completeness.

% Last two columns here are "formal" errors
% l m cosine sine cosine_stddev sine_stddev
potcoffs=nan(nmonths,addmup(Ldata),6);
% Last two columns here are "calibrated" errors
% l m cosine sine
cal_errors=nan(nmonths,addmup(Ldata),4);

% Loop over the months
for index = 1:nmonths 
    % load geopotential coefficients
    fname1=fullfile(ddir1,datanames{index});

    % Open and scan the file (data from both centers is 10 columns)
    fid = fopen(fname1);
    C = textscan(fid,'%s%s%s%s%s%s%s%s%s%s');
    fclose(fid);

    % Only grab the lines for GRCOF2
    Carray = cat(3,C{:});
    I = strmatch('GRCOF2',Carray(:,1,1),'exact');
    Carray = squeeze(Carray(I,1,:));
    
    % Only want columns 2-7, and as format double
    Carray = Carray(:,2:7);
    lmcosi_month=cellfun(@str2num,Carray);
    % This should be addmup(Ldata)
    [m,n] = size(lmcosi_month);
    
    % Change the order of the coefficients so that 
    % order m goes as [0 01 012 0123 ...]
    new_ordering = zeros(m,6);
    revdel=[0 Ldata:-1:0];
    i=1;
    for j=1:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        new_ordering(j,:) = lmcosi_month(k,:);    
        i=i+1;
    end
    lmcosi_month = new_ordering;
    
    % Remove the mean value of the potential i.e. set 0,0 coff = 0
    lmcosi_month(1,3) = 0;
    
    % Make the geopotential relative to the WGS 84 ellipsoid 
    % A bit redundant since we replace 2,0 shortly
    lmcosi_month(4,3) = lmcosi_month(4,3) - j2;
    lmcosi_month(11,3) = lmcosi_month(11,3) - j4;
    
    % Calculate the midpoint of this data span
    monthstart = datenum([str2num(datanames{index}(7:10))...
                        1 str2num(datanames{index}(11:13))]);
    monthend = datenum([str2num(datanames{index}(15:18))...
                        1 str2num(datanames{index}(19:21))]);
    monthmid = (monthstart+monthend)/2;
    thedates(index) = monthmid;
    
    % Now replace the (2,0) coefficient with the SLR value (referenced to
    % WGS84 above).
    % NOTE: This gives a value different than if you used
    % (column3 - column5) from the SLR data file because that data is
    % referenced to an overall mean, not to the WGS 84 ellipsoid.
    where=slrc20(:,1)>monthstart & slrc20(:,1)<monthend;
    if ~any(where)
        % If there is no SLR value within our specific interval, 
        % use the closest value we have
        [~,where]=min(abs(monthmid - slrc20(:,1)));
    end
    % Need to use slrc20(where,2)
    disp(sprintf('C20 was %12.8e now %12.8e',lmcosi_month(4,3),slrc20(where,2)))
    lmcosi_month(4,3)=slrc20(where,2);
    
    % Convert the geopotential coefficients into surface mass density, 
    % if so desired
    if strcmp(units,'SD')
        % Need to make geoid first
        a=fralmanac('a_EGM96','Earth');
        lmcosi_extra=plm2pot([lmcosi_month(:,1:2) lmcosi_month(:,5:6)*a],[],[],[],4);
        lmcosi_month=plm2pot([lmcosi_month(:,1:2) lmcosi_month(:,3:4)*a],[],[],[],4);
        % Add the fornal errors back in columns 5,6
        lmcosi_month=[lmcosi_month lmcosi_extra(:,3:4)];
    end
    
    % Combine into one matrix 
    potcoffs(index,:,:) = lmcosi_month;
    
    % CALIBRATED ERRORS
    fname2=fullfile(ddir1,errornames{index});
    
    % Open and scan the file (data from both Pcenters is 5 columns)
    fid = fopen(fname2);
    E = textscan(fid,'%s%s%s%s%s');
    fclose(fid);

    % Only grab the lines for CALSDV
    Earray = cat(3,E{:});
    I = strmatch('CALSDV',Earray(:,1,1),'exact');
    Earray = squeeze(Earray(I,1,:));
    
    % Only want columns 2-5, and as format double
    Earray = Earray(:,2:5);
    cal_errors_month=cellfun(@str2num,Earray);
    [m,n] = size(cal_errors_month);
    
    % Change the order of the coefficients so that 
    % order m goes as [0 01 012 0123 ...]
    revdel=[0 Ldata:-1:0];
    i=1;
    if Pcenter == 'CSR'
      new_ordering = zeros(m,4);
      demm=dems;
      for j=1:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        new_ordering(j,:) = cal_errors_month(k,:);
        demm(j)=cal_errors_month(k,2);
        i=i+1;
      end
      cal_errors_month = new_ordering;
    elseif strmatch('GSM-2_2006121-2006151_0028_EIGEN_G---_0004.txt',...
                    errornames(index),'exact')
      % for one very odd GFZ file
      new_ordering = zeros(m,4);
      i=4;
      for j=4:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        % This file has only 1 space for the 2,1 coefficients
        if dems(i)==0
          k=k-2;
        else
          k=k-3;
        end
        new_ordering(j-3,:) = cal_errors_month(k,:);
        i=i+1;
      end
      cal_errors_month = new_ordering;
      
    else % for the rest of GFZ, which has slightly less odd formatting
      new_ordering = zeros(m-1,4);
      i=4;
      for j=4:length(dems)
        k = dels(i)+1 + sum( revdel( (1:dems(i) + 1 ) ) );
        % These files have two spaces for the 2,1 coefficients
        if j == 5
          k=k-3;
        else
          k=k-2;
        end
        new_ordering(j-3,:) = cal_errors_month(k,:);
        i=i+1;
      end
      cal_errors_month = new_ordering;
    end
    
    % If from the GFZ data center, add terms for l=0 and 1
    if Pcenter == 'GFZ'
      cal_errors_month = [0 0 0 0; 1 0 0 0; 1 1 0 0; cal_errors_month];
    end
    
    % Replace the C20 error from GRACE with the C20 error from SLR since we
    % used the C20 coefficient from SLR
    disp(sprintf('C20 error was %12.8e now %12.8e',cal_errors_month(4,3),slrc20_error(where)))
    cal_errors_month(4,3)=slrc20_error(where);
    
    % Convert the geopotential error coefficients into surface mass 
    % density, if so desired
    if strcmp(units,'SD')
        % Need to make geoid first
        a=fralmanac('a_EGM96','Earth');
        cal_errors_month=plm2pot([cal_errors_month(:,1:2) cal_errors_month(:,3:4)*a],[],[],[],4);
    end
    
    % Combine into one matrix
    cal_errors(index,:,:) = cal_errors_month;
    
    disp(['Processed month number: ' num2str(index)]);
end

% Save
save(fnpl,'potcoffs','cal_errors','thedates','datanames','errornames');

end % End if we have a save file already

% Collect output
varns={potcoffs,cal_errors,thedates,datanames,errornames};
varargout=varns(1:nargout);

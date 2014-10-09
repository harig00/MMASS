function varargout=raydyntrace(srcloc,recloc,raytype,modelname)
% [rradius,rdelta,ttimeuncor,ellcor,srcloc,recloc]=...
%         RAYDYNTRACE(srcloc,recloc,raytype,modelname)
%
% A simple interface to the FORTRAN program RAYDYNTRACE
%
% INPUT:
% 
% srcloc      [lat lon depth] of source in degrees|km (can be column vectors)
% recloc      [lat lon eleva] of receiver in degrees|km (can be column vectors)
% raytype     'P'  for a P wave [default]
% modelname   'IASP91' for the IASP91 seismic wavespeed model [default]
%
% OUTPUT:
% 
% rradius      Cell array with the radii visited by the ray
% rdelta       Cell array with the epicentral distance progress of the ray
% ttimeuncor   Cell array with the uncorrected travel time of the ray
% ellcor       Cell array with the ellipticity correction of each ray
% srcloc       Source coordinates (from the input)
% recloc       Receiver coordinates (from the input)
%
% EXAMPLE:
%
% [rradius,rdelta,ttimeuncor,ellcor,srcloc,recloc]=raydyntrace;
% circ(6371-[0 660 410]); hold on
% for ind=1:length(rradius)
%    [X,Y]=pol2cart(pi/2-rdelta{ind},rradius{ind});
%    pr(ind)=plot(X,Y);
% end
% axis([-6600 6600 -6600 6600]); axis off
%
% Last modified by fjsimons-at-alum.mit.edu, 10/23/2009

% STEP 1: Retrieve the package from the CIG website link
% http://geoweb.princeton.edu/pub/nolet/Breviary_software.tar
% STEP 2: Compile according to the instructions, set the filename 
% for CRUST20 in the FORTRAN code before you do. See Raydyntrace.manual.
% I used "/usr/local/intel/ifort-10.1/bin/ifort" which worked on uname -a 
% Linux lemaitre 2.6.18-128.2.1.el5 #1 SMP Wed Jul 8 11:54:47 EDT 2009
% x86_64 x86_64 x86_64 GNU/Linux 
% STEP 3: Define paths and environmental variables and run this program

% Define paths, etc - I work from environmental variables
defval('execpath',fullfile(getenv('FFILES'),'BREVIARY','1D-global'))
defval('modelname','IASP91')

% Make up some defaults just in case
defval('srcloc',...
       [79.888  1.856 0.0   ; 79.888    1.856 0.0   ; 79.888   1.856 0.0  ])
defval('recloc',...
       [42.639 74.494 1.645 ; 34.946 -106.457 4.712 ; 25.123 102.740 0.706])
defval('raytype','P')

% Check input parameter formalities
if prod(size(srcloc))~=3 && size(srcloc,2)~=3; balk; end
if prod(size(recloc))~=3 && size(recloc,2)~=3; balk; end
if any(size(srcloc)~=size(recloc)); balk; end

% Produce the appropriate MODEL file
modelfile=fullfile(execpath,'Examples',modelname);
if exist(modelfile,'file')~=2; balk; end

% Produce the appropriate RAY file
switch raytype
 case 'P'
  rayname='Pdef';
  rayfile=fullfile(execpath,'Examples',rayname);
  maxarr=1;
  tdifmax=0;
  kpole=0;
  tableflag=0;
  % The length of the Pdef file in number of lines
  gratreps=5;
 otherwise
  error('This ray type is not supported yet')
end

% Generate some values for all source/receiver combinations
dates=datestr(date,'yyyymmdd');
stations='FJS';

% Produce the appropriate STATION/RECEIVER file
srecname='SRCREC';
srecfile=fullfile(execpath,'Examples',srecname);
fid=fopen(srecfile,'w');
rwfmt='%8s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4s %i\n';
for nrays=1:size(srcloc,1)
  fprintf(fid,rwfmt,...
	  dates,srcloc(nrays,:),recloc(nrays,:),stations,kpole);
end
fclose(fid);

% Produce the appropriate RUN file
runname='RUNFILE';
runfile=fullfile(execpath,'Examples',runname);
fid=fopen(runfile,'w');
fprintf(fid,'%s\n',modelname);
fprintf(fid,'%i %i\n',maxarr,tdifmax);
fprintf(fid,'%s\n',srecname);
fprintf(fid,'%s\n',rayname);
fprintf(fid,'%i\n',tableflag);
fclose(fid);

% Then run the actual code from the directory
olddir=pwd;
cd(fullfile(execpath,'Examples'))
unix(sprintf('%s < %s',...
	     fullfile(execpath,'raydyntrace'),...
	     runfile));

% Now take a look at the output
outname=sprintf('raydyntrace.%s',srecname);
outfile=fullfile(execpath,'Examples',outname);

% First I need to preallocate arrays by knowing how many sets of how many
% nodes are to follow... from the file structure, this should work
[s,t]=unix(sprintf('awk ''$3=="0" && NF==8'' %s',outfile));
% Check that I get what I expect - how many newlines does t have?
nraysexp=sum(abs(t)==10);
difer(nrays-nraysexp,[],[],NaN)
% Peel off the number of ray nodes that we should be expecting
t=parse(t,str2mat(10));
for ind=1:nrays
  rnodesexp{ind}=str2num(deblank(rindeks(parse(t(ind,:),' '),2)));
end

% Initialize cell arrays with this many empties; only do this for the
% ones you really want; add the other ones later
[rradius,rdelta,ttimeuncor,ellcor]=...
    deal(cellnan(nrays,[rnodesexp{:}],ones(1,nrays)));

% Now we are ready for the file read
fid=fopen(outfile);
% First the gratuitous repeats of the inputs
skipf(gratreps,fid);
% Now read, ray per ray, what the file has to offer
for ind=1:nrays
  % Then the one-line-per-segment parameters for each ray
  skipf(gratreps-2,fid)
  % The the actual ray headers - first line
  C=textscan(fgetl(fid),'%s%f%f%f%f%f%f%s');
  % The the actual ray headers - second line
  D=textscan(fgetl(fid),'%f%f%f%f%f%f%f%f');
  % Always one arrival only - check with input
  rmaxarr=D{1}; difer(rmaxarr-maxarr,[],[],NaN)
  % Number of ray nodes to follow - check initial idea
  rnodes=D{2}; difer(rnodes-rnodesexp{ind},[],[],NaN)
  % Travel time in seconds, uncorrected - I WANT THIS
  ttimeuncor{ind}=D{4};
  % Ellipticity correction - I WANT THIS
  ellcor{ind}=D{5};
  % t-star and quality factor
  tstar=D{6}; Q=D{7};
  % Slowness in rad/s
  slownessrads=D{8};

  % Now read in the rays themselves, then parse
  E=textscan(fid,'%f%f%f%f%f%f%f',rnodes);
  % The depth that the ray visits - I WANT THIS
  rradius{ind}=E{1};
  rayanglerad=E{2};
  % Epicentral distance for each segment - I WANT THIS
  rdelta{ind}=E{3};
  velnodekms=E{4};
  oneoverQ=E{5};
  H11seckmmin2=E{6};
  H22seckmmin2=E{7};
  % This is weird here, I can't see it using "cat" but I clearly need it
  skipf(1,fid)
end
fclose(fid);

% Go back to where you came from
cd(olddir)

% Prepare output as desired
varns={rradius,rdelta,ttimeuncor,ellcor,srcloc,recloc};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function balk
error('Input wrong or insufficient')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function skipf(N,fid)
for in=1:N; fgetl(fid); end

% [X1,Y1,Z1]=sph2cart(...
%     srcloc(:,2)*pi/180,srcloc(:,1)*pi/180,6371-srcloc(:,3)); 
% [X2,Y2,Z2]=sph2cart(...
%     recloc(:,2)*pi/180,recloc(:,1)*pi/180,6371-recloc(:,3)); 

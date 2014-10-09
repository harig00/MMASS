function raykernels(dirname,flats)
% RAYKERNELS(dirname,flats)
%
% This is a wrapper routine to another wrapper routine, namely RAYDYNTRACE!
%
% INPUT:
%
% dirname      A directory name string [defaulted]
% flats        Earth flattening factor [defaulted]
% raytype     'P'  for a P wave [default]
% modelname   'IASP91' for the IASP91 seismic wavespeed model [default]
%
% Given a directory with kernel directories (and nothing else), it
% produces the meta files that contain the central ray from source to
% receiver that should correspond to the kernel in the corresponding
% directory. It knows that the Earth is elliptical and uses the same
% Earth flattening factor as in RAYDYNTRACE.f.
%
% Last modified by fjsimons-at-alum.mit.edu, 10/23/2009

% Define default values
defval('dirname','/u/fjsimons/MyPapers/WAVELETS/StephenJudd/kernelsEtc')
% This is the WGS84 value for the flattening of the Earth
defval('flats',1/298.257223563)
defval('raytype','P')
defval('modelname','IASP91')

% The actual conversion factor used -- CLOSE to value in RAYDYNTRACE.f,
% which appears based on an 1864 flattening from Hansen. Now who's old.
flets=(1-flats)^2;

% Get all the directories within this directory
alldirs=ls2cell(dirname);

% Go through them one by one
for index=1:length(alldirs)
  % Reads the meta.data file
  fid=fopen(fullfile(dirname,alldirs{index},'meta.data'));
  % This and ff will need to be adapted according to the meta.data
  skipf(2,fid)
  % This is the true GEOCENTRIC source location string
  srcEtc=textscan(fgetl(fid),'%s%f%f%f');
  % Just making sure
  difer(strcmp(srcEtc{1},'srcEtc')-1,[],[],NaN)
  % This is the true GEOCENTRIC receiver location string
  rcvEtc=textscan(fgetl(fid),'%s%f%f%s%s');
  % Close the file
  fclose(fid);
  % Just making sure
  difer(strcmp(rcvEtc{1},'rcvEtc')-1,[],[],NaN)
  % That's it, now we have what we need
  srclat=srcEtc{2}; % in degrees
  srclon=srcEtc{3}; % in degrees
  srcdep=srcEtc{4}; % depth, in km
  reclat=rcvEtc{2}; % in degrees
  reclon=rcvEtc{3}; % in degrees
  recelv=0; % elevation, in km
  % Now do the reverse mapping to GEOGRAPHIC locations such that, when
  % fed as in input to RAYDYNTRACE, it converts again to GEOCENTRIC, to
  % end up where we started
  srclat=atan(tan(srclat*pi/180)/flets)*180/pi;
  reclat=atan(tan(reclat*pi/180)/flets)*180/pi;
  % Say what epicentral distance this particular example is:
  [epikm,epidelta]=grcdist([srclon srclat],[reclon reclat]);
  disp(sprintf('\n%s Epicentral distance %4.1f\n',alldirs{index},epidelta))
  % And with this, run the program RAYDYNTRACE
  [rradius,rdelta,ttimeuncor,ellcor]=raydyntrace(...
      [srclat srclon srcdep],[reclat reclon recelv]);
  % And now put the output neatly in to a file in the same directory
  fid=fopen(fullfile(dirname,alldirs{index},'ray.data'),'w');
  fprintf(fid,'%f %f\n',[rradius{1} rdelta{1}]');
  fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function skipf(N,fid)
for in=1:N; fgetl(fid); end

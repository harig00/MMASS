function regionalize(regi)
% REGIONALIZE(regi)
%
% Constructs databases of geologic regions from the USGS.
% Later on these get projected and saved with other data.
%
% INPUT: 
% 
% regi    'WORLDGeol1' to do everything, i.e. the master set
%         'AFRGeol1' to do specifically Africa
%         'SAGeol1' to do specifically South America
%         [Other regions to come which aren't from the USGS]
%
% Written by dongwang-at-princeton.edu, 04/16/2010
% Last modified by fjsimons-at-alum.mit.edu, 3/2/2011
%
% Previously known as: make_region_world_1, make_region_Africa_1,
% make_region_SouthAmerica_1.m 

% Where is the data kept or stored?
worldir=fullfile(getenv('IFILES'),'GEOLOGY','WORLD');
% What is the Matlab version of this?
allfname=fullfile(worldir,'WORLDGeol1.mat');
% What is the specific region you want today?
fname=fullfile(worldir,sprintf('%s.mat',regi));

% If the regio doesn't exist you'll need to make it
if exist(fname,'file')~=2
  % If the world doesn't exist you'll need to make that
  if exist(allfname,'file')~=2
    % Load the original shape file data
    geoall=shaperead(fullfile(worldir,'WEP_PRVG'));
    
    % Retrieve useful information - there is much more!
    for ik=1:length(geoall)
      WORLDGeol1(ik).X=geoall(ik).X;
      WORLDGeol1(ik).Y=geoall(ik).Y;
      WORLDGeol1(ik).GEOLPROV=geoall(ik).NAME;
      WORLDGeol1(ik).AREA=geoall(ik).AREA;
    end

    % Data source
    bigweb='http://certmapper.cr.usgs.gov';
    docId='F75EF7E7-963B-413A-9539-BCA88D755B01';
    src1=sprintf('[1] %s/data/we/dds60/spatial/shape/wep_prvg.zip',bigweb);
    src2=sprintf('[2] %s/rooms/utilities/full_metadata.jsp?docId={%s}',bigweb,docId);
    src3=sprintf('[3] http://pubs.usgs.gov/dds/dds-060/');
    source=sprintf('%s \n%s \n%s ',src1,src2,src3);
    
    % Commit to file
    save(allfname,'WORLDGeol1','source')
  else
    load(allfname)
  end
  % Now you make and save the region
  switch regi
   case 'AFRGeol1'
    regionrange=[65 -55 -40 80];
   case 'SAGeol1'
    regionrange=[40 -80 240 360];
    % Change it to [0 360]
    for ik=1:length(WORLDGeol1)
      WORLDGeol1(ik).X(WORLDGeol1(ik).X<0)=...
          WORLDGeol1(ik).X(WORLDGeol1(ik).X<0)+360;
    end
  end
  
  % Now pick out what matters
  jk=0;
  for ik=1:length(WORLDGeol1)
    tempindx=WORLDGeol1(ik).X>=regionrange(3) & ...
             WORLDGeol1(ik).X<=regionrange(4);
    tempindy=WORLDGeol1(ik).Y>=regionrange(2) & ...
             WORLDGeol1(ik).Y<=regionrange(1);
    tempind=tempindx & tempindy;
    
    if ~isempty(find(tempind))
      jk=jk+1;
      eval(sprintf('%s(jk)=WORLDGeol1(ik);',regi));
      eval(sprintf('%s(jk).X(~tempind)=NaN;',regi));
      eval(sprintf('%s(jk).Y(~tempind)=NaN;',regi));
    end
  end
  % SAVING 
  save(fname,regi,'source')
end

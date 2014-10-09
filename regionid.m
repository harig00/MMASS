function varargout=regionid(stro,dfile)
% [af,foundit]=regionid(stro,dfile)
%
% Queries databases of geologic regions from the USGS. 
% 
% INPUT:
%
% stro        Some string, e.g. 'Canar', if empty, lists
% dfile       Some string identifying the database file
%
% EXAMPLE:
%
% [a,b]=regionid('Ocean'); 
% plot([a(:).X],[a(:).Y])
%
% [a,b]=regionid('Valborg','WORLDGeol1.mat');
% 
% Last modified by fjsimons-at-alum.mit.edu, 3/27/2012

% Where is the data kept or stored?
worldir=fullfile(getenv('IFILES'),'GEOLOGY','WORLD');
%worldir=fullfile(getenv('IFILES'),'VENUS');

% What is the Matlab version of this?
defval('dfile',fullfile(worldir,'WORLDGeol1.mat'));

defval('stro',[])
defval('foundit',[])
defval('af',[])

% If the world doesn't exist you'll need to make that
if exist(dfile,'file')~=2
  error('Run REGIONALIZE or at least SHAPEREAD and then resave')
  % Example
  % venus=shaperead('VENUS_nomenclature.shp');
  % Retrieve useful information - there is much more!
  % for ik=1:length(venus)
  %     WORLDGeol1(ik).X=venus(ik).X;
  %     WORLDGeol1(ik).Y=venus(ik).Y;
  %     WORLDGeol1(ik).GEOLPROV=venus(ik).FEATURE;
  %     WORLDGeol1(ik).DIAMATER=venus(ik).DIAMETER;
  %   end    
  %   save(dfile,'WORLDGeol1')
else
  load(dfile)
  a={WORLDGeol1(:).GEOLPROV};
  if ~isempty(stro)
    for in=1:length(a)
      if strfind(a{in},stro); 
        foundit=[foundit in];
      end       
    end
  else
    for in=1:length(a)
      disp(sprintf('%4i %s',in,a{in}))
    end
  end
end

% Return output
if ~isempty(foundit)
  af=WORLDGeol1(foundit);
end
vars={af,foundit};
varargout=vars(1:nargout);

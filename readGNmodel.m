function v=readGNmodel(dep,fname)
% v=READGNMODEL(dep,fname)
%
% Reads a tomographic mantle model saved as an ASCII file by the code
% MAPWAVELET.f made by Guust Nolet and Raffa Montelli and published in
% G-Cubed, 2006. This takes a while, which is not fun.
%
% INPUT:
% 
% dep          Depth requested, in km; if NaN, returns possibilities 
% fname        A full filename to the model in cubed-sphere output
%
% OUTPUT:
%
% v            The NxNx6 values to plot on the cubed sphere, OR:
%              the possible depth values that you may choose
%
% EXAMPLE:
%
% d=readGNmodel(NaN);
% v=readGNmodel(361.200);
% plotoncube(readGNmodel(609.525),'2D')
%
% Last modified by fjsimons-at-alum.mit.edu, 10/25/2010

% Hardwired Earth radius, not ideal
rearthkm=6371;

defval('dep',rearthkm-5761.475)
defval('dep',rearthkm-6009.800)
defval('dep',rearthkm-5716.325)
defval('dep',rearthkm-6167.825)
defval('dep',rearthkm-5942.075)
defval('fname',fullfile(getenv('IFILES'),...
			'EARTHMODELS','MONTELLI','PRI5w'))
defval('xver',0)

% Here you put in what you know ahead about the model - not elegant, but
% sufficient for now
lfin=7;
N=2^lfin+1;

v=[];
fs=fsize(fname);
% disp(sprintf('Opening %s of size %i Mb',suf(fname,'/'),round(fs/1e6))) 
fid=fopen(fname,'r');
% Read header
junk=fgetl(fid);
while 1
  % Load one depth of the model
  onedepth=shiftdim(reshape(fscanf(fid,'%f',N*N*6*4),4,N*N*6),1);
  rads=unique(onedepth(:,3));
  % Check that it is really only one depth
  difer(length(rads)-1,[],1,NaN)

  if isnan(dep)
    % Merely collect the depths assuming the shallowest is first
    v=[v rads];
    if ftell(fid)==fs-1
      v=max(v)-v;
      break
    end
  else
    % Keep reading until you've found the single requested depth
    if abs(dep-(rearthkm-rads))<1e-4
      disp(sprintf('Found requested radius %8.3f or depth %8.3f',...
		   rads,rearthkm-rads))
      % Start building the model
      v=reshape(onedepth(:,4),[N N 6]);
      % Need to figure out what GN convention is with respect to ours
      if xver==1
	lo=reshape(onedepth(:,1),[N N 6]);
	la=reshape(onedepth(:,2),[N N 6]);
	for ind=1:6
	  [xi,eta]=sphere2cube(lo(:,:,ind),la(:,:,ind));
	  [p,pg]=plotonchunk(xi,eta);
	  pause
	end
      end
      % Rearrange in the Judd/Vetter scheme
      v=v(:,:,[5 4 6 2 3 1]);
      fclose(fid);
      break
    end
  end
end


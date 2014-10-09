function varargout=readJRmodel(in,rmav,degres,toplot)
% [cosi,dels,dems]=readJRmodel
% lmcosi=readJRmodel(in,rmav)
% [r,lon,lat]=readJRmodel(in,rmav,degres)
% [r,lon,lat,ch,ph,cb,xl]=readJRmodel(in,rmav,degres,1)
%
% Reads/converts/plots the spherical harmonic coefficients of the seismic
% wavespeed model S40RTS by Jeroen Ritsema.
% 
% INPUT: (if no input, returns the entire model over all splines)
%
% in        The number of the spline "layer" function [default: 3]
% rmav      Remove the degree-zero contribution [default: 0]
% degres    Degree resolution of the expansion [default: 1]
% toplot    Make a Mollweide map [default: 0]
%
% OUTPUT:
%
% cosi,dels, dems  The coefficients for all of the splines/layers, OR
% lmcosi           [degr order cos sin] coefficients for a single layer, OR
% r, lon,lat       The map values and the longitudes and latitudes AND
% ch,ph,cb,xl      Coninent, plate, colorbar, and xlabel handles of the figure
%
% SEE ALSO: PLM2XYZ, PLOTPLM, INTERPJRMODEL, RTSSPLINES, CALCJRMODEL
%
% EXAMPLE:
%
% readJRmodel(3,[],1,1) % should be the same as
% [a,dels,dems]=readJRmodel; plotplm([dels dems a(:,:,3)],[],[],[],1);
% kelicol % Note that this plots a single spline, which isn't yet a "depth"
%
% Last modified by fjsimons-at-alum.mit.edu, 01/19/2011

% Define default file names and directory locations
defval('diro',fullfile(getenv('IFILES'),'EARTHMODELS','RITSEMA'))
defval('fname','S40RTS.sph')

% Open and load the file and read it in
fid=fopen(fullfile(diro,fname));

% Read the first line
fline=fgetl(fid);

% Will parse this line later but now I just "know" it
L=40; n=21;

% Read all data in one go - not sure why it doesn't like %11e4
cofs=reshape(fscanf(fid,'%f'),(L+1)^2,n);
fclose(fid);

% Check that the size is what we expect
difer(prod(size(cofs))-(L+1)^2*n,[],[],NaN)

% Create the vectors with harmonic degrees and orders, etc
[dems,dels,mz,lmcosi,mzi,mzo]=addmon(L);

% Stick the coefficients in at the right place
cosi=repmat(lmcosi(:,3:4),[1 1 n]);
% We're doing this in a roundabout way but the result counts
coso=cosi(:,:,1);
for ind=1:n
  coso(mzo)=cofs(:,ind);
  cosi(:,:,ind)=coso;
end

% Remove the spherical average of the RELATIVE PERTURBATION also
sphav=squeeze(cosi(1,1,:));
defval('rmav',0)
if rmav==1
  cosi(1,1,:)=0;
end

% Need to include the (-1)^m factor also
CSC=(-1).^dems;
% A further adjustment that is needed and was thoroughly verified
dom=1./sqrt(2-(dems==0))/sqrt(4*pi);

if nargin>0
  % Identify the spline/layer to plot
  defval('in',4);
  % Collect the sine and cosine coefficients containing the model
  plotcosi=cosi(:,:,in).*repmat(CSC.*dom,1,2);
  % Change the relative to the PERCENTAGE perturbation
  plotcosi=plotcosi*100;
  
  % Prepare for output 
  varns={[dels dems plotcosi]};

elseif nargin==0 % if nothing else requested return all coefficients
  for in=1:size(cosi,3)
    cosi(:,:,in)=cosi(:,:,in).*repmat(CSC.*dom,1,2)*100;
  end
  
  % Prepare for output 
  varns={cosi,dels,dems};
end

% Output this if there are no further requests 
if nargin>2
  % Perform the expansion
  defval('degres',1)
  [modl,lon,lat]=plm2xyz([dels dems plotcosi],degres);
  
  % Prepare for output 
  varns={modl,lon,lat};
end  

if nargin>3
  defval('toplot',0)
  if toplot==1
    % Do the plotting in Mollweide
    clf ; [d,ch,ph]=plotplm(modl,[],[],1,1);
    % Adjust the colorbar to a symmetric one
    caxis(round(halverange(modl,75)))
    % Make a color scale that is "tomographic"
    kelicol
    % Put a color bar on
    cb=colorbar('hor'); shrink(cb,1.5,1.5);
    longticks(cb,2); axes(cb)
    xl=xlabel('shear velocity variation from 1-D');
    movev(cb,-.05); set(cb,'xaxisl','b')
    
    % Cosmetics
    set(cat(1,ch{:}),'linew',1.5); 
    set(ph,'color','g','linew',1.5)
    fig2print(gcf,'portrait')
    % Prepare for output 
    varns={modl,lon,lat,ch,ph,cb,xl};
  end
end

% Collect output
varargout=varns(1:nargout);
  

function varargout=blocktile(matrix,wlen,operc,indo,xofyof)
% [nwy,nwx]=BLOCKTILE(matrix,wlen,operc)
% tile=BLOCKTILE(matrix,wlen,operc,indo)
% tile=BLOCKTILE(matrix,wlen,operc,indo,[xofset yofset])
% [nwy,nwx]=BLOCKTILE(matrix,wlen,operc,[],[xofset yofset])
%
% Returns (the number of) overlapping squares constituting a data set. 
%
% INPUT:
%
% matrix     A data set
% wlen       Tile length (square tiles)
% operc      Percentage overlap
% indo       Requested running index, down the rows first
% xofyof     Initial offset for first tile
%
% OUTPUT:
%
% nwy        Expected number of tiles in the Y direction
% nwx        Expected number of tiles in the X direction
% tile       The actual tile of the requested index
%
% EXAMPLE 1:
%
% mat=peaks(64);
% subplot(221) ; imagesc(blocktile(mat,30,60,1))
% subplot(222) ; imagesc(blocktile(mat,30,60,3))
% subplot(223) ; imagesc(blocktile(mat,30,60,2))
% subplot(224) ; imagesc(blocktile(mat,30,60,4))
%
% EXAMPLE 2:
%
% blocktile(mat,32,0,3)==blockisolate(mat,[32 32],3)
%
% See also: BLOCKISOLATE
%
% Last modified by fjsimons-at-mit.edu, 28.10.2005

[ny,nx]=size(matrix);

if nargin>4
    [xofset,yofset]=deal(xofyof(1),xofyof(2));
  else
    [xofset,yofset]=deal(0);
end

nx=nx-xofset;
ny=ny-yofset;

% Overlap percentage
operc=floor(operc/100*wlen);
% Number of windows in X
nwix=fix((nx-operc)/(wlen-operc)); 
% Number of windows in Y
nwiy=fix((ny-operc)/(wlen-operc)); 

% If not actual tile requested but the number of them is wanted
if nargin==3 | [nargin==5 & isempty(indo)]
  varargout{1}=nwiy;
  varargout{2}=nwix;
  return
else
  % Which square are we talking about with our running index 'indo'?
  % Transform to matrix coordinates
  xindo=ceil(indo/nwiy);
  yindo=indo-(xindo-1)*nwiy;
  
  colix = 1 + (0:(nwix-1))*(wlen-operc)+xofset;
  coliy = 1 + (0:(nwiy-1))*(wlen-operc)+yofset;
  
  % Calculate and output the tile
  tile=matrix(coliy(yindo):coliy(yindo)+wlen-1,...
      colix(xindo):colix(xindo)+wlen-1);

  varargout{1}=tile;
end




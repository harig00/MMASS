function [bm,xc,yc,CT]=blockmean(mat,side,olap)
% [bm,xc,yc,CT]=BLOCKMEAN(mat,[iside jside],[olapi olapj])
%
% Returns the mean of isideXjside blocks of
% side length, tiled over the matrix.
% If two arguments are given, tiles do not overlap.
% If three arguments are given, overlaps the tiles.
%
% When working with real coordinates, of course these indices
% are pixel-centered. Also returns the center point of the boxes.
% These are the centers for the averaging regions; with the overlap,
% it is not possible to use them as the pixels centers for an image plot.
%
% mat=peaks(64);
% blockmean(mat,[1 1]) is the null operation
% blockmean(mat,size(mat)) returns mean(mat(:))
% blockmean(mat,[4 4],[0 0])==blockmean(mat,[4 4])
%
% mat=rand(120,80);
% for index=1:77
% tile=blocktile(mat,20,50,index);
% difm(index)=abs(mean(tile(:))-indeks(blockmean(mat,[20 20],[10 10]),index));
% end
% disp(sprintf('%8.3e\n',difm))
%
% This is a good function, it seems. Can get rid of first if
% and just make default overlap zero.
%
% Last modified by fjsimons-at-alum.mit.edu, 12/28/2006

[iside,jside]=deal(side(1),side(2));
if nargin==2
  if any(mod(size(mat),side))
    error('Matrix not right size for nonoverlapping tiles')
  end
  % Column space averaging
  ro=1:size(mat,2);
  co=gamini(1:size(mat,2)/jside,jside);
  CT=sparse(ro,co,1);
  post=mat*CT;
  % Row space averaging
  ro=1:size(mat,1);
  co=gamini(1:size(mat,1)/iside,iside);
  CT=sparse(ro,co,1)';
  bm=CT*post;
  bm=bm/prod(side);
  [xc,yc]=deal(NaN);
else
  [olapi,olapj]=deal(olap(1),olap(2));
  [ny,nx]=size(mat);
  nwj=(nx-olapj)/(jside-olapj); % Number of windows in X
  nwi=(ny-olapi)/(iside-olapi); % Number of windows in Y
  if ~all(fix([nwi nwj])==[nwi nwj])
    error('Matrix not right size for overlapping tiles')
  end  
  % Column space averaging
  ro=pauli(1:size(mat,2),jside);
  ro=ro(1:jside-olapj:end,:)';
  ro=ro(:)';
  co=gamini(1:length(ro)/jside,jside);
  CT=sparse(ro,co,1);
  post=mat*CT;
  % Row space averaging
  ro=pauli(1:size(mat,1),iside);
  ro=ro(1:iside-olapi:end,:)';
  ro=ro(:)';
  co=gamini(1:length(ro)/iside,iside);
  CT=sparse(ro,co,1)';
  bm=CT*post;
  bm=bm/prod(side);
  xc=(1+(jside-1)/2):jside-olapj:size(mat,2);
  yc=(1+(iside-1)/2):iside-olapi:size(mat,1);
end


function [mapt,nuh,kurt,sdev]=clutser(matrix,varargin)
% [mapt,nuh,kurt,sdev]=CLUTSER(matrix)
% [mapt,nuh,kurt,sdev]=CLUTSER(matrix,ml1,ml2,ml3)
%
% After this, plot
% bar(nuh(2,:),nuh(1,:))
%
% Returns mapping matrix that maps the new matrix
% into the old one, trying to make a flatter histogram
% of its entries - compare nuh and kurt.
% ml1, ml2 and ml3 are empirical values - 
% multiplying the mean of the input matrix; if that product is
% higher than what you get over a square then do not subdivide futher.
%
% See also NEQPLOT.
%
% Last modified by fjsimons-at-alum.mit.edu, 09/11/2007

% Largest cell is 8x8, smallest cell is 1x1 in this case

% Do this for every block of 8X8 elements
if any(mod(size(matrix),8))
  error('Matrix size not a multiple of 8')
end

mm=mean(matrix(:));
mapt=repmat(NaN,size(matrix,1),size(matrix,2));
maxi=0;

wb=waitbar(0,'Progress...');

for index=1:prod(size(matrix)/8)
  waitbar(index/prod(size(matrix)/8),wb)
  [mat,i,j]=blockisolate(matrix,[8 8],index);
  % Big enough to subdivide further?
  bm8=blockmean(mat,[8 8]); % This is of course the mean of the isolated matrix
  bm4=blockmean(mat,[4 4]);
  bm2=blockmean(mat,[2 2]);
  % Criterium: if the function over the blocks has a higher mean
  % than the overal mean of the matrix, then accept subdivision.
  if nargin==1
    ml1=7; ml2=3; ml3=1;
  else
    [ml1,ml2,ml3]=deal(varargin{1},varargin{2},varargin{3});
  end  
  crit22=ml1*mm;  crit44=ml2*mm; crit88=ml3*mm;
  bm8=gamini2(bm8,[8 8])>=crit88; % If 1 accept 4x4; If 0 stay 8x8
  bm4=gamini2(bm4,[4 4])>=crit44; % If 1 accept 2x2; If 0 stay 4x4
  bm2=gamini2(bm2,[2 2])>=crit22; % If 1 accept 1x1; If 0 stay 2x2

  % Now take care criteria that are compatible:
  % division needs to be cascading
  big=bm2.*bm4.*bm8+bm4.*bm8+bm8;
  
  sm1=size(mat,1);
  sm2=size(mat,2);

  % Takes care of column space
  ind1=repmat([1:4 ; sm1+[1:4]],1,sm2/2)+gamini2([0:sm2/2-1]*2*sm1,[2 4]);
  
  % Now take care of the row space
  ind1=repmat(ind1,1,sm1/4)+gamini2([0:sm1/4-1]*4,[2 2*sm2]);
  
  % Map ready for cumsum
  % Make 0 into 1/64, 1 into 1/16, 2 into 1/4 and 3 into 1
  big=1./4.^(3-big); 
  big(ind1(:))=ceil(cumsum(big(ind1(:))))+maxi;
  maxi=max(big(:));
% For testing:  maxi=0;
  mapt(i,j)=big;
end
close(wb)

% Now look at distribution; compare hist(matrix) with hist(matrix(mapt))
nudens=sparse(1,mapt(:),matrix(:));
[nuh,binm]=hist(nudens(~~nudens),25);
sdev=full(std(nudens(~~nudens))); 
kurt=kurtosis(nudens(~~nudens));
nuh=[nuh ; binm];

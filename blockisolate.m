function [bloke,i,j]=blockisolate(matrix,nmsize,indo)
% [bloke,i,j]=BLOCKISOLATE(matrix,[n m],indo)
%
% Isolates the 'indo'-th block matrix of size
% nXm out of a 'matrix' 
%
% EXAMPLE:
%
% mat=peaks(64);
% subplot(221) ; imagesc(blockisolate(mat,[32 32],1))
% subplot(222) ; imagesc(blockisolate(mat,[32 32],3))
% subplot(223) ; imagesc(blockisolate(mat,[32 32],2))
% subplot(224) ; imagesc(blockisolate(mat,[32 32],4))
%
% See also BLOCKTILE, BLOCKISOLATE2
%
% Written by fjsimons-at-mit.edu, Oct 23th 2001

[n,m]=deal(nmsize(1),nmsize(2));

if any(mod(size(matrix),nmsize))
  error('Matrix size not a multiple of n or m')
end

mulm=size(matrix,1)/m;
muln=size(matrix,2)/n;

jindo=ceil(indo/mulm);
iindo=indo-(jindo-1)*mulm;

i=(iindo-1)*m+1:iindo*m;
j=(jindo-1)*n+1:jindo*n;

bloke=matrix(i,j);


function [th,p,scl,avHes,Fisher,truecov,nh]=osrzero(fid,np,npp)
% [th,p,scl,avHes,Fisher,truecov,nh]=OSRZERO(fid,np,npp)
%
% Reads a THZRO file as we have come to use it 
%
% INPUT:
%
% fid        The file id (from FOPEN)
% np         The number of parameters in the vector
% npp        The number of nonzero elements in the matrices 
%
% OUTPUT:
% 
% All the good stuff
%
% SEE ALSO: 
%
% OSWZEROE, OSWZEROB, OSLOAD
%
% Last modified by fjsimons-at-alum.mit.edu, 10/06/2014

fgetl(fid);
% The truth
th=fscanf(fid,'%f',np)';

% The other parameters of the experiment, see SIMULOS and MLEOS
fgetl(fid); fgetl(fid);

% If this involves gravity - or not
if np>=5
  p=fscanf(fid,'%f',9)';
else
  p=fscanf(fid,'%f',5)';
end

% The scale used in the optimization
% Only once for some legacy code 
fgetl(fid); fgetl(fid);
scl=fscanf(fid,'%f',np)';

% This is the theoretical covariance of the estimate
fgetl(fid); fgetl(fid);
truecov=fscanf(fid,'%f',npp)';

% This is the right scaled Fisher matrix from which the above derives 
fgetl(fid); fgetl(fid);
Fisher=fscanf(fid,'%f',npp)';

% And the average Hessian could be close to the Fisher if you're lucky
% Don't necessarily look at THIS partial average of the Hessians
% through the iterations as it's just of the last few iterations that
% add cumulatively to THINI and THHAT. We are getting this again from
% the full file diagnos. So if we have interrupted a sequence of
% simulations we need run one more simulation to close out this
% file properly, which we do by setting N=0 in MLEROS etc.
fgetl(fid); 
% Pick out the number that got reported
t=fgetl(fid); nh=str2num(t([abs(t)<58 & abs(t)>47]));
avHes=fscanf(fid,'%f',npp)';

% Now reassemble the covariance
fullcov=zeros(np,np);
fullcov(nonzeros(triu(reshape(1:np^2,np,np)')'))=truecov;
truecov=[tril(fullcov)'+tril(fullcov)-diag(diag(fullcov))];

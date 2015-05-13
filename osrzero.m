function [th,p,scl,avHes,Fisher,covF,nh]=osrzero(fid,np)
% [th,p,scl,scl,avHes,Fisher,covF,nh]=OSRZERO(fid,np)
%
% Reads a THZRO file as produced by the olhede? suite of
% programs following Simons & Olhede (2013).
%
% INPUT:
%
% fid        The file id (from FOPEN)
% np         The number of parameters in the vector
%
% OUTPUT:
% 
% th         The true parameter vector
% p          The parameter structure of the simulation
% scl        The scaling factors
% avHes      The average Hessian matrix
% Fisher     The unblurred Fisher matrix
% covF       The predicted covariance matrix, based on Fisher
% nh         The number of simulations yielding avHes
%
% SEE ALSO: 
%
% OSWZEROE, OSWZEROB, OSLOAD
%
% Last modified by fjsimons-at-alum.mit.edu, 02/09/2015

% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% Start the read
fgetl(fid);
% The unscaled truth
th=fscanf(fid,'%f',np)';

% The other parameters of the experiment, see SIMULOS and MLEOS
fgetl(fid); fgetl(fid);

% If this involves gravity - or not
if np>=5
  fields={'DEL','g','z2','dydx','NyNx','blurs','kiso'};
  p=fscanf(fid,'%f',10)';
  valjus={[p(1:2)] p(3) p(4) [p(5:6)] p(7:8) p(9) p(10)};
else
  fields={               'dydx','NyNx','blurs','kiso'};
  p=fscanf(fid,'%f',6)';
  valjus={[p(1:2)] [p(3:4)] p(5) p(6)};
end
% Structurize those values
p=cell2struct(valjus,fields,2);

% Here you need to read the optimization options until you're done
fgetl(fid); fgetl(fid);
g=':';
while ~isempty(strfind(g,':'))
  g=fgetl(fid);
end

% Here you need to read the optimization bounds until you're done
fgetl(fid);
g=':';
while ~isempty(strfind(g,':'))
  g=fgetl(fid);
end

% The scale used for the Fisher matrix
scl=fscanf(fid,'%f',np)';

% This is the theoretical covariance of the estimate
fgetl(fid); fgetl(fid);
covF=trilosi(fscanf(fid,'%f',npp));

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


function x=preconD6(x,precon,tipe,cofs)
% x=preconD6(x,precon,tipe,cofs)
%
% Preconditions a three-dimensional array for D6 interval
% transforms. Interval wavelets are merely orthogonal; the
% preconditioning takes care of preserving the moments.
%
% INPUT:
%
% x         The three-dimensional array, dimensions must be powers of two
% precon    Array identifying preconditioning, e.g. [1|0 1|0] to
%           precondition or not in the first two (angular) dimensions
% tipe      'forward'|'inverse'|'transpose'|'inversetranspose'
% cofs      The filter coefficients coming out of D6BOXCOF
%
% OUTPUT:
%
% x         The preconditioned output array
%
% Last modified by fjsimons-at-alum.mit.edu, 10/13/2010

defval('xver',0)

% Figure out dimensions
nall=size(x);

if length(nall)<2; nall(2)=1; end
if length(nall)<3; nall(3)=1; end
if length(precon)<2 ; precon(2)=0; end
if length(precon)<3 ; precon(3)=0; end

% In which non-singleton dimensions do we precondition?
dims=find([precon~=0] & [nall~=1]);

% Determine the coefficients - let's fix the peculiar transfos in the
% coef file later
switch tipe 
 case 'forward'
  tops=cofs.LF';
  bots=flipud(fliplr(cofs.RF));
 case 'transpose'
  tops=cofs.LF;
  bots=flipud(flipud(cofs.RF)');
 case 'inversetranspose'
  tops=cofs.LI;
  bots=flipud(flipud(cofs.RI)');
 case 'inverse'
  tops=cofs.LI';
  bots=flipud(fliplr(cofs.RI));
end

% Then actually do it
for dim=dims
  % Isolate the sets of planes in the right dimension
  xleft= [x(dindeks(          1,dim,nall))'; ...
          x(dindeks(          2,dim,nall))'; ...
          x(dindeks(          3,dim,nall))'];
  xright=[x(dindeks(nall(dim)-2,dim,nall))'; ...
	  x(dindeks(nall(dim)-1,dim,nall))'; ...
    	  x(dindeks(nall(dim)-0,dim,nall))'];

  % Perform the actual preconditioning in the right dimension
  stuff=tops*xleft;
  for i=1:length(cofs.H0)/2
    x(dindeks(          i,dim,nall))=stuff(i,:);
  end

  stuff=bots*xright;
  for i=1:length(cofs.H0)/2
    x(dindeks(nall(dim)-(length(cofs.H0)/2-i),dim,nall))=stuff(i,:);
  end

  if xver==1
    disp(sprintf('%s preconditioned in dimension %i',tipe,dim))
  end
end

% If nothing was done, keep the original



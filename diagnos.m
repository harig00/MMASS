function [thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covh]=diagnos(fname,ddir,np)
% [thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covh]=DIAGNOS(fname,ddir,np)
%
% Reads in a single file with diagnostics from MLEOS, MLEROS0, MLEROS.
%
% INPUT:
%
% ddir    A directory name, e.g. OLHEDE?/Simulations-lemaitre
% fname   A file name, e.g. diagn_01-Jul-2014
% np      The number of parameters that are being solved for
%
% OUTPUT:
%
% See MLEOS, MLEROS0
%
% thhat    The maximum-likelihood estimate of the vector with elements:
%          [       s2 nu rho], in "nothing", see SIMULOSL
%          [D f2 r s2 nu rho], in Nm, and "nothing", see SIMULROS
%          [D f2   s2 nu rho], in Nm, and "nothing", see SIMULOS, SIMULROS0
% thini    The starting guess used in the optimization
% tseiter  Time taken, error code, number of iterations
% scl      The scaling applied as part of the optimization procedure
% L        The likelihood of this solution
% gam      The score, or first derivative of the likelihood
% hes      The Hessian, second derivative of the likelihood
% optis    The first-order optimality condition returned by FMINUNC
% momx     The various moments of the quadratic piece of the likelihood
% covh     The asymptotic covariances derived from the Hessian
%
% SEE ALSO:
%
% OSOPEN, OSLOAD (with which it needs to match!)
%
% Last modified by fjsimons-at-alum.mit.edu, 02/11/2015

defval('ddir','/u/fjsimons/PROGRAMS/MFILES/olhede4')
defval('fname','mleosl_diagn_02-Oct-2014')

% The number of parameters to solve for; standard is 5
defval('np',5)
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% The number of sample variances that will be read - reset below
defval('nvar',2)

% Get a rough estimate of the number of estimates from the size
% You want this number to err on the side of being too large!
if np==3
  nsize=450; nvar=1;
elseif np==5
  nsize=560;
elseif np==6
  nsize=721;
end
ndim=ceil(fsize(fullfile(ddir,fname))/nsize);

% Initialize
tseiter=nan(ndim,3);
L=nan(ndim,1);
optis=nan(ndim,1);
momx=nan(ndim,3);
% Make room for the data sample variance(s)
thhat=deal(nan(ndim,np+nvar));
[thini,gam]=deal(nan(ndim,np));
hes=nan(ndim,npp);
covh=nan(ndim,npp);

% Rarely, in SPMD mode does the file get written too quickly and does a
% confusion between labs happen - look into the file and fix easily
% Read the contents
fid=fopen(fullfile(ddir,fname),'r');
for index=1:ndim
  try
    % The initial guess
    thini(index,:)=fscanf(fid,'%e',np);
  catch
    % Quit if you're done
    index=index-1;
    break
  end
  % The estimates, and the scaled spatial sample variance(s)
  thhat(index,:)=fscanf(fid,'%e',np+nvar);
  % Three diagnostics (time taken, exitflag, number of iterations)
  tseiter(index,:)=fscanf(fid,'%i',3); 
  try
    % The likelihood
    L(index)=fscanf(fid,'%e',1);L(index)
  catch
    % Tell if you're wrong
    L(index-1)
    sprintf('%12.8g',thhat(index-1,3))
    break
  end
  % The first-order optimality criterion
  optis(index)=fscanf(fid,'%e',1);
  % The two moments of Xk, and the magic parameter
  momx(index,:)=fscanf(fid,'%e',3);
  % The scalings
  scl(index,:)=fscanf(fid,'%e',np);
  % The scores
  gam(index,:)=fscanf(fid,'%e',np);
  % The scaled Hessian elements at the solution
  hes(index,:)=fscanf(fid,'%f',npp);
  % The unscaled covariance from the Hessian at the solution
  covh(index,:)=fscanf(fid,'%f',npp);
end
fclose(fid);

% Finetune the preallocation over time
disp(sprintf('Expected length %i, true length %i',size(thhat,1),index))

% Trim the possibly wrongly preallocated arrays
thhat=thhat(1:index,:);
thini=thini(1:index,:);
tseiter=tseiter(1:index,:);
scl=scl(1:index,:);
L=L(1:index,:);
gam=gam(1:index,:);
hes=hes(1:index,:);
covh=covh(1:index,:);
optis=optis(1:index,:);
momx=momx(1:index,:);

% Put out
varns={thhat,thini,tseiter,scl,L,gam,hes,optis,momx,covh};
varargout=varns(1:nargout);

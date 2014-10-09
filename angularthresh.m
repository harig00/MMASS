function [vwtinv,vw,vwt,recer,kilme]=...
    angularthresh(v,tperc,J,iface,wav,precon,xver)
% [vwtinv,vw,vwt,recer,kilme]=...
%   ANGULARTHRESH(v,tperc,J,iface,wav,precon,xver)
%
% Hard thresholding of a model in cubed-sphere coordinates
% 
% INPUT:
%
% v         The NxNx6 input model
% tperc     Truncation level as percentile [defaulted]
% J         Maximum scale in both directions
% iface     Numbered face if you don't want to do them all [0]
% wav       'D2', 'D4' or 'D6' identifying the basis
% precon    [1|0 1|0] toggles preconditioning in either dimension 
% xver      1 Provides optional truth and fact checking
%           0 Doesn't [default]
% 
% OUTPUT:
% 
% vwtinv    The NxNx6 output model after thresholded reconstruction
% vw        The original wavelet and scaling coefficients
% vwt       The thresholded wavelet and scaling coefficients
% recer     An array with the following diagnostics:
%           [1] The reconstruction error l2-NORM, in percent
%           [2] The number of zeroed wavelet coefficients, in percent 
%           [3] The actual truncation value being used
%           [4] The truncation level as a percentile
%           [5] The l1 norm of the thresholded coefficients, in percent
%           [6] The total variation of the reconstruction, in percent
% kilme     The logicals identifying where the wavelets are killed
%
% Last modified by fjsimons-at-alum.mit.edu, 1/21/2011

defval('iface',0)
defval('tperc',85)
defval('precon',[1 1])
defval('xver',0);

if iface>0
  v=v(:,:,iface);
end

switch wav
 case 'D2'
  % The 2-tap Daubechies wavelet transform
  vw=angularD2WT(v,[J J],'forward',1);
 case 'D4'
  % The 4-tap Daubechies wavelet transform
  vw=angularD4WT(v,[J J],precon,'forward',1);
 case 'D6'
  % The 6-tap Daubechies wavelet transform
  vw=angularD6WT(v,[J J],precon,'forward',1);
end

% Test Ignace's code as far as D4 and D6 are concerned
whereitsat=fullfile(getenv('MFILES'),'ignaceloris');
if xver==1 && ~strcmp(wav,'D2')
  disp(sprintf('\n%s\n',...
	       'Extra verification using Ignace Loris'' code'))
  switch wav
   case 'D4'
    addpath(whereitsat)
    difer([vw-angularD4MRWT(v,[J J],precon,'forward')]/prod(size(v)),8); 
    rmpath(whereitsat)
    % Make sure to check my own inverse transform
    difer([v-angularD4WT(vw,[J J],precon,'inverse',1)]/prod(size(v)),8);
   case 'D6'
    addpath(whereitsat)
    difer([vw-angularD6MRWT(v,[J J],precon,'forward')]/prod(size(v)),8);
    rmpath(whereitsat)
    % Make sure to check my own inverse transform
    difer([v-angularD6WT(vw,[J J],precon,'inverse',1)]/prod(size(v)),8);
  end
  disp(sprintf('\n%s\n',...
	       'Done with extra verification using Loris'' code'))
end

% Explicit and absolute zeroing of the MAGNITUDE
trunx=prctile(abs(vw(:)),tperc);
kilme=abs(vw)<trunx;

% With this truncation, how well do you approximate Africa?
vwt=vw;
vwt(kilme)=0;

% In other words, we still need the inversion
switch wav
 case 'D2'
  vwtinv=angularD2WT(vwt,[J J],'inverse',1);
 case 'D4'
  vwtinv=angularD4WT(vwt,[J J],precon,'inverse',1);
 case 'D6'
  vwtinv=angularD6WT(vwt,[J J],precon,'inverse',1);
end

% Now compute the Total Variation
% Take a quick look:
% plotoncube(vwtinv,'3D',[],[],[],[],[],[],0)
tv0=0; tv=0;
for i=1:size(v,3)
  [Dx,Dy]=gradient(v(:,:,i));
  % Take note, this is like ell_1 of the gradient!
  tv0=tv0+sum(sum(sqrt([Dx.^2+Dy.^2])));
  [Dx,Dy]=gradient(vwtinv(:,:,i));
  tv =tv +sum(sum(sqrt([Dx.^2+Dy.^2])));
end

% The relative reconstruction error, etc, see the help
recer=[norm(vwtinv(:)-v(:))/norm(v(:))*100 ...
       sum(kilme(:))/prod(size(kilme))*100 ...
       trunx ...
       tperc ...
       norm(vwt(:),1)/norm(vw(:),1)*100 ...
       tv/tv0*100];

function varargout=chunkwt(x,precon,tipe,levels,bases)
% f=CHUNKWT(x,precon,tipe,[n1 n2 n3],{'b1' 'b2' 'b3'})
%
% Performs a wavelet transform on ONE OF six cubed-sphere chunks. With
% proper accounting for the edges of the box. It doesn't see like the
% multiresolution properties are well handled here. See other transforms.
%
% INPUT:
%
% x             The three-dimensional input array 
% precon        Array identifying three preconditioning flags, 1 or 0
% tipe          'forward'|'inverse'|'transpose'|'inversetranspose'
% [n1 n2 n3]    The number of levels in each direction (default: [3 3 3])
% 'b1' etc      Strings identifying the particular wavelet basis in the
%                two azimuthal and one radial dimension, e.g.
%               {'D4' 'D4' 'D2'} where D is for Daubechies and 6, 4, 2 are
%               allowed for the number of TAPS, 2 being the Haar
%               transform, and the first ones needing to be identical
%
% OUTPUT:
%
% f             The size(x) output array
% wc            The cell array with the indices to the wavelet coefficients
% sc            The cell array with the indices to the scaling coefficients
% 
% SEE ALSO:
%
% ANGULARD4WT, ANGULARD6WT, D4BOXSTEP, D4BOXSTEPI, D6BOXSTEP, D6BOXSTEPI
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2010

% Set some default values
defval('x',rand(128,128,37))
defval('precon',[0 0 0])
defval('tipe','forward')
defval('levels',[3 3 3])
defval('bases',{'D4' 'D4' 'D2'})

% Assign the levels
n1=levels(1);
n2=levels(2);
n3=levels(3);

% Assign the bases
b1=bases{1};
b2=bases{2};
b3=bases{3};

% Repeat the defaults just in case
defval('b1','D4')
defval('b2','D4')
defval('b3','D2')

% Now perform the transform in the AZIMUTHAL dimensions
if strcmp(b1,'D4') && strcmp(b2,'D4')
  % It's what we've always done:
  % D4 in the first two dimensions...
  f=angularD4WT(x,levels(1:2),precon(1:2),tipe);
  % ...and D2 aka Haar in the third, radial, dimension
  f=radialHaarWT(f,levels(3),precon(3),tipe);
elseif strcmp(b1,'D6') && strcmp(b2,'D6')
  % It's the new think
  % D6 in the first two dimensions...
  f=angularD6WT(x,levels(1:2),precon(1:2),tipe);
else 
  error('Specify valid options (D4|D6) for the azimuthal transforms')
end
% And perform the transform in the RADIAL dimensions
if strcmp(b3,'D2') 
  % ...and D2 aka Haar in the third, radial, dimension
  f=radialHaarWT(f,levels(3),precon(3),tipe);
end

% Prepare the output
varns={f};
varargout=varns(1:nargout);


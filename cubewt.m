function F=cubewt(X,precon,tipe,levels,bases)
% F=CUBEWT(X,precon,tipe,[n1 n2 n3],{'b1' 'b2' 'b3'})
%
% Performs a wavelet transform on A SET OF six cubed-sphere chunks. With
% proper accounting for the edges of the box.
%
% INPUT:
%
% X           The input structure with six three-dimensional arrays, as 
%                X.xplus, X.zminus, X.yplus, X.xminus, X.zplus, X.yminus 
% precon      Array identifying three preconditioning flags, 1 or 0
% tipe        'forward'|'inverse'|'transpose'|'inversetranspose'
% [n1 n2 n3]  The number of levels in each direction (default: [3 3 3])
% 'b1' etc    Strings identifying the particular wavelet basis, e.g.
%             {'D6' 'D4' 'D2'} where D is for Daubechies and 6, 4, 2 are
%             for the number of TAPS, thus 2 being the Haar transform
% 
% OUTPUT:
%
% F             The output array with the wavelet and scaling coefficients
% 
% SEE ALSO:
%
% CHUNKWT, RANDCUBE, NANCUBE, PLOTONCUBE3
%
% EXAMPLE:
%
% cubewt('demo1')
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2010

if ~isstr(X)
  % Supply the defaults
  defval('X',randcube)
  defval('precon',[0 0 0])
  defval('tipe','forward')
  defval('levels',[4 4 4])
  defval('bases',{'D4' 'D4' 'D2'})

  % Initialize the output array
  fnX=fieldnames(X);
  % This assumes all cubes have identical dimensions, of course
  szX=size(X.(fnX{1}));
  F=nancube(szX(1),szX(2),szX(3),fnX);

  % Do the transform
  for index=1:length(fnX)
    disp(sprintf('Working on face %s',fnX{index}))
    F.(fnX{index})=chunkwt(X.(fnX{index}),precon,tipe,levels,bases);
  end
elseif strcmp(X,'demo1')
  F=randcube;
  FF=cubewt(F,[],'forward',[4 4 4]);
  FI=cubewt(FF,[],'inverse',[4 4 4]);
  difer(F-FI)
  % Plot some with PLOTONCUBE2
end

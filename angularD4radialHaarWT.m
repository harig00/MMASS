function varargout=angularD4radialHaarWT(x,levels,precon,meth,tipe,mr)
% f=angularD4radialHaarWT(x,levels,precon,meth,tipe,mr)
%
% Performs D4 wavelet transform in the angular direction (i.e. in the
% first and second index), and Haar wavelet transform in the radial
% direction (i.e. in the third index).
%
% INPUT:
%
% x            The three-dimensional array
% levels       [n1,n2,n3],
%              the number of levels in the three directions
%              n1=0, n2=0, n3=0 is the identity if precon=[0 0 0]
% precon       Array of length 3 identifying preconditioning  
%      [1 1 1] Precondition (NOT an orthogonal transform) so a sequence of
%              ones and a linearly increasing sequence are mapped to a
%              sequence of zeroes in the wavelet bands 
%      [0 0 0] Don't precondition (YES, this is an orthogonal transform)
%              so a sequence of ones and a linearly increasing sequence are
%              not mapped to all zeroes at the edges in the wavelet bands,
%              i.e. the first and last two coefficients are not zero
% meth  'pixel'     Angular reduction by numerical integration
%       'wavelets'  Angular reduction by wavelet transformation
% tipe         'forward'|'inverse'|'transpose'|'inversetranspose', always of
%              the same meth in both indices
% mr         0 The "original" improper, dimensionally sequential transform
%            1 The "proper" multiresolution transform [default]
%
% OUTPUT:
% 
% f            The size(x) matrix with the transform coefficients; scaling
%              coefficients in front, followed by wavelet coefficients. 
%
% EXAMPLE:
%
% angularD4radialHaarWT('demo1')
% angularD4radialHaarWT('demo2')
%
% Last modified by fjsimons-at-alum.mit.edu, Nov 7, 2011
% Last modified by yanhuay-at-princeton.edu, Nov 16, 2011

if ~isstr(x)
  % Get the levels of the decomposition
  defval('levels',[3 3 3])
  defval('meth', 'pixel')
  defval('tipe','forward')
  defval('mr',1)
  defval('precon',[1 1 1])
  
  % Determine the original size of the input
  [neta,nxi,nrad]=size(x);
  
  % Radial reduction from neta*nxi*129 to neta*nxi*37
  if nrad~=37
    x=radialpart129to37(x);
  end

  % Angular reduction from 513*513*nrad to 128*128*nrad
  if neta~=128 || nxi~=128
    x=angularpart513to128(x,meth);
  end
  
  % Haar wavelet transform in radial direction 
  levelHaar=levels(3);
  preconHaar=precon(3);
  x=radialHaarWT(x,levelHaar,preconHaar,tipe);
  
  % D4 wavelet in angular direction
  levelD4=levels(1:2);
  % If the first reduction is by wavelets, no more need for
  % preconditioning, see IL document (2/7/2009) page 7
  preconD4=[1 1]*strcmp(meth,'pixel');
  % Do the angular transform
  x=angularD4WT(x,levelD4,preconD4,tipe,mr);
  
  % Prepare the output
  varns={x};
  varargout=varns(1:nargout);
  
elseif strcmp(x,'demo1')
  disp('Check whether inverse is inverse of forward:')
  num1=513;
  num2=513;
  numr=129;
  for lev=3:3
    x=randn(num1,num2,numr);
    fp=angularD4radialHaarWT(x,[lev lev lev],[],'pixel','forward');
    xp=angularD4radialHaarWT(fp,[lev lev lev],[],'pixel','inverse');

    fw=angularD4radialHaarWT(x,[lev lev lev],[],'wavelets','forward');
    xw=angularD4radialHaarWT(fw,[lev lev lev],[],'wavelets','inverse');
    %Test the error 
    % need to consider whether we could inverse reduction xp/xw back to size of x  
    if numr~=37
      x37=radialpart129to37(x);
    else
      x37=x;
    end
    if num1~=128 || num2~=128
      xsmallp=angularpart513to128(x37,'pixel');
      xsmallw=angularpart513to128(x37,'wavelets');
    else
      xsmallp=x37;
      xsmallw=x37;
    end
    
    errp=max(max(max(abs(xp-xsmallp))));
    disp(['Levels=' num2str(lev) ' The maximal error for pixel method is ' num2str(errp)]);
    if(errp<1e-6)
      disp('Wavelet transform can be completely constructed')
    else
      disp('Wavelet transform can not be completely constructed')
    end
    disp(' ')
    
    errw=max(max(max(abs(xw-xsmallw))));
    disp(['Levels=' num2str(lev) ' The maximal error for wavelet method is ' num2str(errw)]);
    if(errw<1e-6)
      disp('Wavelet transform can be completely constructed')
    else
      disp('Wavelet transform can not be completely constructed')
    end
    disp(' ')
  end
  
elseif strcmp(x,'demo2') 
  disp('Check the inner product:')
  num1=513;
  num2=513;
  numr=129;
  for lev=3
    % Construct a fine kernel
    kernel=ones(num1,num2,numr);
    % By pixel method
    kwav1=angularD4radialHaarWT(kernel,[lev lev lev],[],'pixel','inversetranspose');
    model=ones(128,128,37); 
    modelwav1=angularD4radialHaarWT(model,[lev lev lev],[],'pixel','forward');
    p1=dotproduct3d(modelwav1,kwav1);
    
    % By wavelet
      % reduction
    kwav2=angularD4radialHaarWT(kernel,[lev lev lev],[],'wavelets','inversetranspose');
    model=ones(512,512,37); 
    model=angularD4WT(model,[2 2],[1 1],'forward');
    model=model(1:128,1:128,:);
    modelwav2=angularD4radialHaarWT(model,[lev lev lev],[],'wavelets','forward');
    p2=dotproduct3d(modelwav2,kwav2); 
    errpw=max(max(max(abs((p1-p2)/p1))));
    disp(['kernel=constant. Constant model. Integral err is:' num2str(errpw)]);
    disp(' ')   
  end
end


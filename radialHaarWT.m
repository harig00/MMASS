function varargout=radialHaarWT(x,nz,precon,tipe)
% varargout=RADIALHAARWT(x,nz,precon,tipe)
%
% Perform Haar wavelet transform in the third index of a
% three-dimensional array of size MxNx37 (hardwired) up to nz scales.
%
% INPUT:
%
% x         The three-dimensional array of size MxNx37
% nz        The number of levels in the third direction; 
%           nz=0 takes care of thin layers (NOT identity)
%           nz must be smaller than 5 for 37 layers
% precon    1 Precondition (NOT an orthogonal tranformation) so a
%             sequence of ones in the radial direction is mapped to a
%             sequence of zeroes in the wavelet bands and a sequence of
%             ones in the low pass band under the forward transform 
%           0 Don't precondition (YES, this is an orthogonal transform)
%             so a sequence of ones is not mapped to all zeroes in
%             wavelet bands (thinner layers)
% tipe      'forward'|'inverse'|'transpose'|'inversetranspose'
% 
% OUTPUT:
% 
% f          The size(x) matrix with the transform coefficients; scaling
%            coefficients in front, followed by wavelet coefficients. 
%
% EXAMPLE:
%
% radialHaarWT('demoX'), where X=1,2,3,4,5
%
% Written by Ignace Loris (igloris@vub.ac.be) on 22/06/2009
% Last modified by fjsimons-at-alum.mit.edu, 02/18/2010
% Modified by yanhuay@princeton.edu 11/16/2011

if ~isstr(x)
  defval('nz',3)
  defval('tipe','forward')
  defval('precon',1)

  % Error check on hardwired dimensions
  if size(x,3)~=37 
    error('Input should have size 37 in third variable!');
  end
  if nz>5 
    error('Maximum number of Haar wavelet scales is 5!');
  end

  % Number of thick layers
  numb=32;

  % Preconditioning
  a=1/sqrt(2);
  if precon
    % disp('Preconditioning! No longer an orthogonal transform')
    % Need to pre-multiply by these numbers 
    % Wavelet coefficients of constant will then be zero
    % Note that in the following, a for-loop is actually faster
    % Compare this with the layer thicknesses in this special
    % construction knowing that the first index is the deepest level
    prefactor=repmat([a a ones(1,23) a a 1 1 a a 1 1 1 a 1/2 1/2],...
		     [size(x,2) 1 size(x,1)]);
    % Special provision if x is effecively unidimensionsl
    if size(x,1)~=1 && size(x,2)~=1
      prefactor=shiftdim(prefactor,2);
    else
      prefactor=reshape(prefactor,size(x));
    end
    switch tipe
     case 'forward'
      x=x.*prefactor;
     case 'inversetranspose'
      x=x./prefactor;
     case {'transpose','inverse'}
     otherwise
      error('specify valid type')
    end
  end

  % The following are orthogonal transformations
  switch tipe
   case {'forward','inversetranspose'}
    % Initialize
    f=x;
    % Take care of thin layers (work from the bottom up!)
    % Put averages in correct place
    f(:,:, 1)=a*(x(:,:, 1)+x(:,:, 2));
    f(:,:,25)=a*(x(:,:,26)+x(:,:,27));
    f(:,:,28)=a*(x(:,:,30)+x(:,:,31));
    f(:,:,32)=(x(:,:,36)+x(:,:,37))/2+a*x(:,:,35);
    
    % Put differences at the back
    f(:,:,33)=a*(x(:,:, 1)-x(:,:, 2));
    f(:,:,34)=a*(x(:,:,26)-x(:,:,27));
    f(:,:,35)=a*(x(:,:,30)-x(:,:,31));
    f(:,:,36)=a*x(:,:,35)-(x(:,:,36)+x(:,:,37))/2;
    f(:,:,37)=a*(x(:,:,36)-x(:,:,37));

    % Move normal layers in correct position
    f(:,:,[2:24 26:27 29:31])=x(:,:,[3:25 28:29 32:34]);

    % Now you have 32 equally spaced ~90 km thick layers
    
    % Now do nz of iterations
    for level=1:nz
      k=numb/2^level;
      x=f;
      %lowpass
      f(:,:,1:k)    =a*(x(:,:,1:2:2*k-1)+x(:,:,2:2:2*k));
      %highpass
      f(:,:,k+1:2*k)=a*(x(:,:,1:2:2*k-1)-x(:,:,2:2:2*k));
    end
    
   case {'inverse','transpose'}
    f=x;
    % First undo nz of iterations
    for level=nz:-1:1
      x=f;
      k=numb/2^level;
      f(:,:,1:2:2*k-1)=a*(x(:,:,1:k)+x(:,:,k+1:2*k));
      f(:,:,2:2:2*k)  =a*(x(:,:,1:k)-x(:,:,k+1:2*k));
    end
    x=f;
    % Put the thick layers in the right place
    f(:,:,[3:25 26 28 29 32:34])=x(:,:,[2:24 25:27 29:31]);
    % Now reconstruct the thin layers
    f(:,:, 1)=a*(x(:,:, 1)+x(:,:,33));
    f(:,:, 2)=a*(x(:,:, 1)-x(:,:,33));
    f(:,:,26)=a*(x(:,:,25)+x(:,:,34));
    f(:,:,27)=a*(x(:,:,25)-x(:,:,34));
    f(:,:,30)=a*(x(:,:,28)+x(:,:,35));
    f(:,:,31)=a*(x(:,:,28)-x(:,:,35));
    f(:,:,35)=a*(x(:,:,32)+x(:,:,36));
    x(:,:,36)=(-x(:,:,36)+a*f(:,:,35))/a;
    f(:,:,36)=a*(x(:,:,36)+x(:,:,37));
    f(:,:,37)=a*(x(:,:,36)-x(:,:,37));
   otherwise
    error('specify valid type')
  end

  % preconditioning
  if precon
    % Need to post-multiply 
    switch tipe 
     case 'inverse'
      f=f./prefactor;
     case 'transpose'
      f=f.*prefactor;
     case {'forward','inversetranspose'}
     otherwise
      error('specify valid type')
    end
  end

  % Optional output
  varns={f};
  varargout=varns(1:nargout);
  
elseif strcmp(x,'demo1')
  %check whether data can be constructed completely
  disp('Check whether inverse is inverse of forward:')
  precon=0; % 0 or 1 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  num1=128;
  num2=128;
  %fixed
  numr=37;
  for lev=0:4
    x=randn(num1,num2,numr);
    f=radialHaarWT(x,lev,precon,'forward');
    y=radialHaarWT(f,lev,precon,'inverse');
    err=max(max(max(x-y)));
    disp(['Levels=' num2str(lev) ' The maximal error is ' num2str(err)]);
    if(err<1e-8)
      disp('Haar wavelet can be completely constructed')
    else
      disp('Haar wavelet can not be completely constructed')
    end
    %disp(' ')
  end
  precon=1; % 0 or 1 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  num1=128;
  num2=128;
  numr=37;%fixed
  for lev=0:4
    x=randn(num1,num2,numr);
    f=radialHaarWT(x,lev,precon,'forward');
    y=radialHaarWT(f,lev,precon,'inverse');
    err=max(max(max(x-y)));
    disp(['Levels=' num2str(lev) ' The maximal error is ' num2str(err)]);
    if(err<1e-8)
      disp('Haar wavelet can be completely constructed')
    else
      disp('Haar wavelet can not be completely constructed')
    end
    %disp(' ')
  end
  
elseif strcmp(x,'demo2')
  %Check orthogonality of Haar wavelet transform, i.e. whether norm can
  %be preserved
  disp('Check whether norm can be preserved:')
  disp(' ')
  precon=0;
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  for lev=0:4
    x=rand(1,1,37);
    disp(['lev=' num2str(lev)])
    % test orthogonality of forward transform, no preconditioning
    w=radialHaarWT(x,lev,precon,'forward');
    ratio=norm3d(x)/norm3d(w);
    disp(['Check orthogonality of forward (no precond). Ratio=' num2str(ratio)])
    % test orthogonality of inverse transform, no preconditioning
    w=rand(1,1,37);
    x=radialHaarWT(w,lev,precon,'inverse');
    ratio=norm3d(x)/norm3d(w);
    disp(['Check orthogonality of inverse (no precond). Ratio=' num2str(ratio)])
    %disp(' ')
  end
  
  precon=1;
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  for lev=0:4
    x=rand(1,1,37);
    disp(['lev=' num2str(lev)])
    % test orthogonality of forward transform, no preconditioning
    w=radialHaarWT(x,lev,precon,'forward');
    ratio=norm3d(x)/norm3d(w);
    disp(['Check orthogonality of forward (no precond). Ratio=' num2str(ratio)])
    % test orthogonality of inverse transform, no preconditioning
    w=rand(1,1,37);
    x=radialHaarWT(w,lev,precon,'inverse');
    ratio=norm3d(x)/norm3d(w);
    disp(['Check orthogonality of inverse (no precond). Ratio=' num2str(ratio)])
    %disp(' ');
  end
  
elseif strcmp(x,'demo3')
  %Check whether inner product can be preserved
  %check whether inner product <x, y>=<Wx, W^(-1,t) y>
  disp('Check inner product:')
  disp('Check whether <x, y>=<Wx, W^(-1,t) y>:')
  precon=1; % 0 or 1 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  num1=128;
  num2=128;
  numr=37;%fixed
  for lev=0:4
    vec1=randn(num1,num2,numr);
    wav2=randn(num1,num2,numr);
    wav1=radialHaarWT(vec1,lev,precon,'forward');
    vec2=radialHaarWT(wav2,lev,precon,'inversetranspose');
    p1=dotproduct3d(vec1,wav2);
    p2=dotproduct3d(wav1,vec2);
    err=(p1-p2)/p1;
    disp(['Inner product error (forward/inversetranspose), levels=' num2str(lev) ': ' num2str(err)]);
    %disp(' ')
  end
  
  %check whether inner product <Wx, y>=<x, W^t y>
  disp('Check whether <Wx, y>=<x, W^t y>:')
  precon=1; % 0 or 1 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  num1=128;
  num2=128;
  numr=37;%fixed
  for lev=0:4
    vec1=randn(num1,num2,numr);
    wav2=randn(num1,num2,numr);
    wav1=radialHaarWT(vec1,lev,precon,'forward');
    vec2=radialHaarWT(wav2,lev,precon,'transpose');
    p1=dotproduct3d(wav1,wav2);
    p2=dotproduct3d(vec1,vec2);
    err=(p1-p2)/p1;
    disp(['Inner product error (forward/transpose), levels=' num2str(lev) ': ' num2str(err)]);
    %disp(' ')
  end

  %check whether inner product <Wx, y>=<x, W^(t,-1) y>
  disp('Check whether <W^(-1)x, y>=<x, W^(t,-1) y>:')
  precon=1; % 0 or 1 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  num1=128;
  num2=128;
  numr=37;%fixed
  for lev=0:4
    vec1=randn(num1,num2,numr);
    wav2=randn(num1,num2,numr);
    wav1=radialHaarWT(vec1,lev,precon,'inverse');
    vec2=radialHaarWT(wav2,lev,precon,'inversetranspose');
    p1=dotproduct3d(wav1,wav2);
    p2=dotproduct3d(vec1,vec2);
    err=(p1-p2)/p1;
    disp(['Inner product error (inverse/inversetranspose), levels=' num2str(lev) ': ' num2str(err)]);
    %disp(' ')
  end     

  %check whether inner product <Wx, y>=<x, W^-1 y>
  disp('Check whether <Wx, y>=<x, W^(-1) y>:')
  precon=1; % precon=0, when precon=1, error is large? 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  num1=128;
  num2=128;
  numr=37;%fixed
  for lev=0:4
    vec1=randn(num1,num2,numr);
    wav2=randn(num1,num2,numr);
    wav1=radialHaarWT(vec1,lev,precon,'forward');
    vec2=radialHaarWT(wav2,lev,precon,'inverse');
    p1=dotproduct3d(wav1,wav2);
    p2=dotproduct3d(vec1,vec2);
    err=(p1-p2)/p1;
    disp(['Inner product error (forward/inverse), levels=' num2str(lev) ': ' num2str(err)]);
    %disp(' ')
  end
  
elseif strcmp(x,'demo4')
  %with preconditioning, wavelet transform can map a constant vector to
  %all zero in the wavlet bands(ref to Ignace page2)
  disp('Check whether constant vector can be mapped to all zeros in the wavelet bands:')
  disp(' ')
  precon=0; 
  if(~precon)
    disp('Without preconditioning:')
    disp('The coeff. in coarse band are not all equal, and those in finer bands are not all zero:')
    disp(' ')
  else
    disp('With preconditioning:')  
    disp('The coeff. in coarse band are all equal, and those in finer bands are all zero:')
    disp(' ')
  end
  num1=1;
  num2=1;
  numr=37;%fixed
  vec=ones(num1,num2,numr); % constant vector
  for lev=4:4
    wav=radialHaarWT(vec,lev,precon,'forward');
    disp(['levels= ', num2str(lev) ' wavelet coefficients: '])
    squeeze(wav)'
    disp(' ')
  end
  
  precon=1; % precon=0 or 1
  if(~precon)
    disp('Without preconditioning:')
    disp('The coeff. in coarse band are not all equal, and those in finer bands are not all zero:')
    disp(' ')
  else
    disp('With preconditioning:')  
    disp('The coeff. in coarse band are all equal, and those in finer bands are all zero:')
    disp(' ')
  end
  num1=1;
  num2=1;
  numr=37;%fixed
  vec=ones(num1,num2,numr); % constant vector
  for lev=4:4
    wav=radialHaarWT(vec,lev,precon,'forward');
    disp(['levels= ' num2str(lev) ', wavelet coefficients: '])
    squeeze(wav)'
    disp(' ')
  end 
elseif strcmp(x,'demo5')
  % now we test that the integral of a product of 
  % 2 functions is equal to the inner product of wavelet coefficients.
  % function 1 (the kernel) starts on 129 radial grid points
  % function 2 (the model) lives on 37 layers from the start
  precon=0; % 0 or 1 
  if(~precon)
    disp('Without preconditioning')
  else
    disp('With preconditioning')   
  end
  
  for lev=0:4
    disp(['lev=' num2str(lev)])
    kernel=ones(1,1,129);
    %incorporate the jacobian and reduces to 37 layers
    kernel37=radialpart129to37(kernel);
    kernelwav=radialHaarWT(kernel37,lev,precon,'inversetranspose');
    model=ones(1,1,37);
    modelwav=radialHaarWT(model,lev,precon,'forward');
    %dot product along the radial index=integral over r
    int2=dotproduct3d(kernelwav,modelwav);
    % The Earth's total radius
    R_E=6371;
    % The Earth's core radius (the core-mantle-boundary)
    R_CMB=3481.4;
    int1=(R_E^3-R_CMB^3)/3; %int1=volume/(4*pi)
    errint=abs(int1-int2)/int1;
    disp(['Integral of product of constant functions. Error=' num2str(errint)])
    disp(' ')
  end
  % now we test that the integral of a product of 
  % 2 functions is equal to the inner product of wavelet coefficients.
  % function 1 (the kernel) starts on 129 radial grid points
  % function 2 (the model) lives on 37 layers from the start
  for lev=0:4
    disp(['lev=' num2str(lev)])
    kernel=ones(1,1,129);
    % The location of the interfaces themselves
    r=linspace(R_CMB,R_E,129);
    Jdr=r.^2;
    % Some trial function that cancels the Jacobian
    kernel(1,1,:)=(1./Jdr);%/delta_R;
    kernel37=radialpart129to37(kernel);
    kernelwav=radialHaarWT(kernel37,lev,precon,'inversetranspose');
    model=ones(1,1,37);
    modelwav=radialHaarWT(model,lev,precon,'forward');
    %dot product along the radial index=integral over r
    int2=dotproduct3d(kernelwav,modelwav);
    %int1=thickness from CMB to Earth surface
    int1=R_E-R_CMB;
    errint=abs(int1-int2)/int1;
    disp(['Integral of product 1/J and constant. Error=' num2str(errint)])
    disp(' ')
  end
end


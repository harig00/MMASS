function xp=angularpart513to128(x,tipe)
% xp=ANGULARPART513TO128(x,tipe)
%
% INPUT:
%
% x        The three-dimensional array of size 513x513xQ
% tipe     'pixel' by numerical integration [fast]
%          'wavelets' by wavelet transformation [slow]
%
% OUTPUT:
%
% xp       The three-dimensional output array of size 128x128xQ
%
% EXAMPLE:
%
% angularpart513to128('demo1')
%
% Written by Ignace Loris (igloris@vub.ac.be) on 23/06/2009
% Last modified by fjsimons-at-alum.mit.edu, 11/07/2011

if ~isstr(x)
  % Error check on hardwired dimensions
  [neta,nxi,nrad]=size(x);
  
  if neta~=513 || nxi~=513
    error('Input should have size 513 in first and second dimension!')
  end

  % Prepare output array
  xp=zeros(128,128,nrad);

  % Produce the Jacobian and include the sampling interval
  [J,coordd,dxi,deta]=cubejac(neta,nxi);
  J=J.*dxi.*deta;

  % For every depth interval
  for i=1:nrad
    % Multiply by Jacobian and sampling interval
    x(:,:,i)=x(:,:,i).*J;
    switch tipe
     case 'pixel'
      % Use Simpson's rule, with weights [1 4 2 4 1]/3;
      % First in xi direction
      y=zeros(128,513);
      y=(1/3)*x(1:4:509,:,i)+(4/3)*x(2:4:510,:,i)+(2/3)*x(3:4:511,:,i)+...
	(4/3)*x(4:4:512,:,i)+(1/3)*x(5:4:513,:,i);
      % Then in eta direction
      xp(:,:,i)=(1/3)*y(:,1:4:509)+(4/3)*y(:,2:4:510)+(2/3)*y(:,3:4:511)+...
		(4/3)*y(:,4:4:512)+(1/3)*y(:,5:4:513);
      % Here is how FJS would have done it... Numerical Recipes 4.1.4 vs 4.1.13
      % tst=randi(128); difer([simpson(1:5,x(1:5,tst,i))-y(1,tst)]/y(1,tst))
     case 'wavelets'
      % Change the dimensions from 513x513xQ to 512x512xQ by averaging
      temp=[x(1:512,1:512,i)+x(2:513,1:512,i)+x(1:512,2:513,i)+x(2:513,2:513,i)]/4;
      % Compute two levels of the preconditioned inverse transpose wavelet
      % transform
      tempw=angularD4WT(temp,[2 2],[1 1],'inversetranspose');
      % Throw away the wavelet subbands
      xp(:,:,i)=tempw(1:128,1:128);
     otherwise
      error('The allowable options are ''pixel'' or ''wavelets''')
    end
  end
elseif strcmp(x,'demo1')
  % Surface area test
  kernel=ones(513,513,1);
  f=angularpart513to128(kernel,'pixel');
  err=abs(sum(sum(f))-4*pi/6)/(4*pi/6);
  disp(['kernel=constant. Model=constant. Integral err is:' num2str(err)])

  % test integral of sin(xi)^2*J(xi,eta)
  x=ones(1,513)'*(sin(linspace(-pi/4,pi/4,513)).^2);
  kernel(:,:,1)=x;
  f=angularpart513to128(kernel,'pixel');
  int=sum(sum(f));
  realint=-sqrt(3)+2*pi/3;
  err=abs(int-realint)/realint;
  disp(['kernel=sin(xi)^2. Model=constant. Integral err is:' num2str(err)])


  % test integral of xi^2*J(xi,eta)
  x=ones(1,513)'*(linspace(-pi/4,pi/4,513).^2);
  kernel(:,:,1)=x;
  f=angularpart513to128(kernel,'pixel');
  int=sum(sum(f));
  realint=0.40882059352406189314;
  err=abs(int-realint)/realint;
  disp(['kernel=xi^2. Model=constant. Integral err is:' num2str(err)])


  % test integral of sin(eta)^2*cos(xi)^2*J(xi,eta)
  x=exp(linspace(-pi/4,pi/4,513))'*ones(1,513);
  kernel(:,:,1)=x;
  f=angularpart513to128(kernel,'pixel');
  int=sum(sum(f));
  realint=2.3050494965745932318;
  err=abs(int-realint)/realint;
  disp(['kernel=sin(eta)^2*cos(xi)^2. Model=constant. Integral err is:' num2str(err)])
elseif strcmp(x,'demo2')
  % Surface area test
  % precompute the coarse wavelet representation of the constant model:
  % constant model
  modelfull=ones(512,512,1);                        
  % reduction
  temp=angularD4WT(modelfull,[2 2],[1 1],'forward');
  modelsmall=temp(1:128,1:128,:);

  % test integral of constant kernel * constant model
  kernel=ones(513,513,1);
  f=angularpart513to128(kernel,'wavelets');
  int=sum(sum(f.*modelsmall));
  err=abs(int-4*pi/6)/(4*pi/6);
  disp(['kernel=constant. Constant model. Integral err is:' num2str(err)])

  % test integral of sin(xi)^2*J(xi,eta)
  x=ones(1,513)'*(sin(linspace(-pi/4,pi/4,513)).^2);
  kernel(:,:,1)=x;
  f=angularpart513to128(kernel,'wavelets');
  int=sum(sum(f.*modelsmall));
  realint=-sqrt(3)+2*pi/3;
  err=abs(int-realint)/realint;
  disp(['kernel=sin(xi)^2. Constant model. Integral err is:' num2str(err)])


  % test integral of xi^2*J(xi,eta)
  x=ones(1,513)'*(linspace(-pi/4,pi/4,513).^2);
  kernel(:,:,1)=x;
  f=angularpart513to128(kernel,'wavelets');
  int=sum(sum(f.*modelsmall));
  realint=0.40882059352406189314;
  err=abs(int-realint)/realint;
  disp(['kernel=xi^2. Constant model. Integral err is:' num2str(err)])



  % test integral of sin(eta)^2*cos(xi)^2*J(xi,eta)
  x=exp(linspace(-pi/4,pi/4,513))'*ones(1,513);
  kernel(:,:,1)=x;
  f=angularpart513to128(kernel,'wavelets');
  int=sum(sum(f.*modelsmall));
  realint=2.3050494965745932318;
  err=abs(int-realint)/realint;
  disp(['kernel=sin(eta)^2*cos(xi)^2. Constant model. Integral err is:' num2str(err)])

end

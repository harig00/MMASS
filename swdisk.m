function varargout=swdisk(m,N,K,L,x,method,scalem)
% [E,V,Nm,c,x,E2]=SWDISK(m,N,K,L,x,method,scalem)
%
% Calculates the radial part of the Slepian functions concentrated to a
% circular domain in Cartesian space. For comparison, can solve the
% integral equation directly, and we can also compare with fully
% two-dimensional methods using LOCALIZATION2D. Note that increasing L is
% almost always a good idea. Remember the distinction between the overall
% Shannon number N and the fixed-order Shannon number Nm, and see NSUBM.
%
% INPUT:
%
% m        The fixed order of the eigenvalue problem
% N        The Shannon number of the full concentration problem [default: 20]
% K        The number of functions (best first) that you want [default: N]
% L        The maximum degree in the Jacobi/Bessel expansion [default: 2N]
%          or else the maximum polynomial degree for which GL will be exact
% x        The abscissas required [default: 0-1 for 'DV' and 0-2 for 'SE']
%          If NaN only computes the eigenvalues
% method   'DV' by the method of De Villiers
%          'SE' by Slepian "extension" [default]
%          'GL' by direct Gauss-Legendre integration
% scalem   1 Scales the solution for weightless orthogonality [default]
%          0 Scales the solution for orthogonality with weight x
%
% OUTPUT:
%
% E        The orthogonal radial Slepian eigenfunctions for the disk,
%          where 2\pi\int_0^1 E_i(x) E_j(x)  \,dx=\delta_{ij}V_i OR
%                2\pi\int_0^1 E_i(x) E_j(x) x\,dx=\delta_{ij}V_in
% V        The concentration eigenvalues
% N        The partial Shannon number for this order
% c        The concentration factor c=2*sqrt(N)
% x        The abscissas at which the functions are evaluated
% E2       The orthonormal radial Slepian eigenfunctions for the disk,
%          where \int_0^1 E_i(x) E_j(x)  \,dx=\delta_{ij} OR
%                \int_0^1 E_i(x) E_j(x) x\,dx=\delta_{ij}
%
% SEE ALSO:
%
% JACOBI, LOCALIZATION2D, KERNELC2D, DEVILLIERS
%
% EXAMPLE:
%
% swdisk('demoX') where X=1->10
%
% Last modified by dongwang-at-princeton.edu, 08/06/2008
% Last modified by fjsimons-at-alum.mit.edu, 03/04/2013

% Supply defaults
defval('m',3)
defval('N',3)
defval('K',ceil(N))
defval('method','SE')
defval('scalem',1)

% Some heuristic defaults - make identical in DEVILLIERS
defval('L',min(max([2*K ceil(2*N) 2*m 10]),...
	       84+192*sum([method=='GL'])))

% Time the code
t0=clock;

if ~isstr(m)  
  switch method
   case 'DV'
    disp('Jacobi expansion method on the interval')
    
    % Get the Jacobi expansion coefficients
    [D,V,~,c,~,~,D2,V2]=devilliers(m,N,K,L);
    % Quote the maximum degree of the Jacobi expansion
    disp(sprintf('\nJacobi expansion up to degree %i\n',L))

    % What is the interval over which we are looking?
    defval('x',linspace(0,1,200));

    if ~isnan(x)
      if any(x<0)
	error('Arguments must be zero or positive')
      end
      
      if ~~sum(x(:)>1 | x(:)<0); 
	warning('Results will be inaccurate outside [0 1]'); 
      end
      
      % Calculate the required set of scaled Jacobi polynomials...
      % ... for the normalization...
      [w,xGL,NGL]=gausslegendrecof(L,[],[0 1]);
      [PGL,PGL2]=getscaledjacobis(xGL,m,L);
      
      % Calculate the required set of scaled Jacobi polynomials...
      % ... for the actual points requested
      [P,P2]=getscaledjacobis(x,m,L);
      
      % Calculate the whole set functions as a matrix with row dimension the
      % evaluation points and column dimension the eigenfunction rank.
      EGL=PGL'*D;
      E=P'*D;

      % This from the symmetrized version
      EGL2=PGL2'*D2;
      ED2=P2'*D2;
      % End of trial with the D2. 

      % Calculate the normalization on the unit interval
      orv2 =2*pi*EGL'*diag(w)*EGL;
      orv22=2*pi*EGL2'*diag(w)*EGL2;

      disp(sprintf('\nGauss-Legendre integration, Nystrom method with %i nodes',...
		   NGL))
      disp(sprintf('Normalization: Mean absolute error %8.3e',...
		   mean(mean(abs(orv2-diag(diag(orv2)))))))
      disp(sprintf('Normalization: Mean absolute error %8.3e',...
		   mean(mean(abs(orv22-diag(diag(orv22)))))))
      % Sometimes you won't want this, and undoing might be hard when the
      % eigenvalues are very small, they may be inaccurate
      E2=sqrt(2*pi)*[diag(1./sqrt(diag(orv2)))*E']';
      % Apply and return the normalized result
      E=[diag(sqrt(V(:))./sqrt(diag(orv2)))*E']';
      %  This wasn't there before, an oversight?
      EGL=[diag(sqrt(V(:))./sqrt(diag(orv2)))*EGL']';
      
      % Now for the symmetrized version
      E22=sqrt(2*pi)*[diag(1./sqrt(diag(orv22)))*ED2']';
      ED2=[diag(sqrt(V2(:))./sqrt(diag(orv22)))*ED2']';
      
      % They may appear in a different order
      difer(sort(E2(:))-sort(E22(:)),[],[],NaN)
    else
      E=NaN; E2=NaN;
    end

   case 'SE'
    % Slepian's method which is valid everywhere
    disp(sprintf('\nBessel expansion method, uniformly valid\n'))

    % Get the Jacobi expansion coefficients
    [D,V,g,c,~,~,D2,V2,g2]=devilliers(m,N,K,L);

    % What is the interval over which we are looking?
    defval('x',linspace(0,1,200));

    if ~isnan(x)
      if any(x<0)
	error('Arguments must be zero or positive')
      end
      
      % Calculate the required set of Bessel functions...
      % ... for the normalization...
      [w,xGL,NGL]=gausslegendrecof(L,[],[0 1]);
      [JGL,JGL2]=getscaledbessels(xGL,m,L,c);

      % ... for the actual points requested
      [J,J2]=getscaledbessels(x,m,L,c);
      
      % Calculate the whole set functions as a matrix with row dimension the
      % evaluation points and column dimension the eigenfunction rank...
      % ... at the integration points
      EGL=JGL'*D./repmat(g,length(xGL),1);
      EGL2=JGL2'*D2./repmat(g2,length(xGL),1);
      
      % ... at the evaluation points
      E=J'*D./repmat(g,length(x),1);
      ED2=J2'*D2./repmat(g2,length(x),1);
      % Prepare the normalization on the unit interval
      orv2=2*pi*EGL'*diag(w)*EGL;
      orv22=2*pi*EGL2'*diag(w)*EGL2;

      disp(sprintf('\nGauss-Legendre normalization, Nystrom method with %i nodes',...
		   NGL))
      disp(sprintf('Normalization: Mean absolute error %8.3e',...
		   mean(mean(abs(orv2-diag(diag(orv2)))))))
      disp(sprintf('Normalization: Mean absolute error %8.3e',...
		   mean(mean(abs(orv22-diag(diag(orv22)))))))
      % Sometimes you won't want this, and undoing might be hard
      E2=sqrt(2*pi)*[diag(1./sqrt(diag(orv2)))*E']';

      % Apply and return the normalized result
      E=[diag(sqrt(V(:))./sqrt(diag(orv2)))*E']';
      EGL=[diag(sqrt(V(:))./sqrt(diag(orv2)))*EGL']';
      
      % Now for the symmetrized version
      E22=sqrt(2*pi)*[diag(1./sqrt(diag(orv22)))*ED2']';
      ED2=[diag(sqrt(V2(:))./sqrt(diag(orv22)))*ED2']';
      
      % They may appear in a different order and their sign may differ
      difer(sort(abs(E2(:)))-sort(abs(E22(:))),[],[],NaN)
    else
      E=NaN; E2=NaN;
    end
       
   case 'GL'
    if isnan(x)
      error('If you only want eigenvalues select ''DV'' or ''SE''')
    end
    % By the Nystrom method on the Bessel kernel
    c=2*sqrt(N);

    % Get the weights and nodes for Gauss-Legendre integration
    [w,xGL,NGL]=gausslegendrecof(L,[],[0 1]);
    
    disp(sprintf('\nGauss-Legendre integration, Nystrom method with %i nodes',...
		 NGL))

    % What is the interval over which we are looking?
    defval('x',linspace(0,1,200));

    if any(x<0)
      error('Arguments must be zero or positive')
    end
    
    switch scalem
     case 1 % More like Slepian, Devilliers and Shkolnisky
      % Evaluate the kernel at the GL points and the requested points
      cxxGL=c*xGL(:)*xGL(:)';
      cxxint=c*x(:)*xGL(:)';
      % Form the actual integration kernel
      KGL=sqrt(cxxGL).*besselj(m,cxxGL);
      Kint=sqrt(cxxint).*besselj(m,cxxint);
     case 0 % More like Simons, Dahlen and Wieczorek
      % This leads to the difference at the origin, right?
      % Ahem - new stuff to be closer to our own without needing
      % 1/sqrt(x) - no "symmetrizing" as Slepian (1964) eq. (19)
      cxxGL=c*xGL(:)*xGL(:)';
      cxxint=c*x(:)*xGL(:)';
      KGL=sqrt(c)*besselj(m,cxxGL).*repmat(xGL(:)',length(xGL),1);
      Kint=sqrt(c)*besselj(m,cxxint).*repmat(xGL(:)',length(x),1);
    end
    
    % Get eigenfunctions of the "symmetrized" kernel - note that Matlab
    % may still consider them numerically asymmetric. And note that the
    % number of nodes appears to be a sensitive parameter.
    if K>=NGL-1 & K<=NGL
      [EGL,g]=eig(diag(sqrt(w))*KGL*diag(sqrt(w)));
    elseif K<NGL-1
      OPTS.disp=0;
      [EGL,g]=eigs(diag(sqrt(w))*KGL*diag(sqrt(w)),K,'LM',OPTS);
    elseif K>NGL
      error(sprintf('You must increase polynomial degree to at least %i',2*K-2))
    end
    
    % Need to unscale E
    EGL=EGL./repmat(sqrt(w(:)),1,size(EGL,2));

    % Now get the actual eigenvalues of the concentration problem
    V=diag(g).^2*c;
    
    % Sometimes you won't want this, and undoing might be hard
    EGL2=EGL;
    % Make sure the inner product over the domain is as wished    
    EGL=EGL*sqrt(diag(V))/sqrt(2*pi);
    
    % Order actual eigenvalues and eigenfunctions downgoing
    [V,i]=sort(V,'descend'); V=V(:)';
    EGL=EGL(:,i); EGL2=EGL2(:,i);

    % But usually we "interpolate" them using the same kernel at the
    % desired values.
    E=(Kint*diag(w))*EGL*inv(g);
    E2=(Kint*diag(w))*EGL2*inv(g);
    
    % Should build in some kind of complex/negative protection
    if any(imag(V)>eps)
      error('Complex eigenvalues - try increasing L');
    else
      V=real(V); E=real(E); E2=real(E2);
    end
    
    % Check that the eigenvalues are properly contained
    if ~isnan(nanmean([V(V>1) V(V<0)]));
      V1=V-1; V1=V1(V1>0); V0=V+1; V0=V0(V0<1);
      if nanmean([V1 V0])>1000*eps
	error('Eigenvalues exceeding 0 to 1 - definitely increase L');
      end
    end

    if scalem==1 || ~strcmp(method,'GL')
      % Check orthogonality over the interval 
      orv2=2*pi*EGL'*diag(w)*EGL;
      checksit1=mean(mean(abs(orv2-diag(V))));
      % Check orthonormality over the interval
      orv3=EGL2'*diag(w)*EGL2;
      checksit2=mean(mean(abs(orv3-eye(size(orv3)))));
    elseif scalem==0 && strcmp(method,'GL')
      orv2=2*pi*EGL'*diag(w)*[EGL.*repmat(xGL,1,size(EGL,2))];
      % Apply and return the normalized result
      EGL=[diag(sqrt(V(:))./sqrt(diag(orv2)))*EGL']';
      EGL2=[diag(sqrt(V(:))./sqrt(diag(orv2)))*EGL2']';
      E=[diag(sqrt(V(:))./sqrt(diag(orv2)))*E']';
      orv2=2*pi*EGL'*diag(w)*[EGL.*repmat(xGL,1,size(EGL,2))];
      checksit1=mean(mean(abs(orv2-diag(V))));
      % Check orthonormality over the interval
      orv3=EGL2'*diag(w)*[EGL2.*repmat(xGL,1,size(EGL2,2))];
      checksit2=mean(mean(abs(orv3-eye(size(orv3)))));
    else
      error('You have run out of options, buddy')
    end
    disp(sprintf('Orthogonality  of E : Mean absolute error %8.3e',...
		 checksit1))
    disp(sprintf('Orthonormality of E2: Mean absolute error %8.3e',...
		 checksit2))
    difer(checksit1,[],[],NaN)
    difer(checksit2,[],[],NaN)
   otherwise
    error('Specify valid method')
  end
  
  % End of the various calculation methods

  if ~strcmp(method,'GL') & ~isnan(x)
    % Check orthogonality over the interval 
    orv2=2*pi*EGL'*diag(w)*EGL;
    % Sort of equal to 2*pi*nansum(E(x<=1,3).^2).*indeks(diff(x))
    % Sort of equal to 2*pi*trapeze(x(x(:)<=1),[E(x<=1,3)].^2)
    checksit1=mean(mean(abs(orv2-diag(V))));
    disp(sprintf('Orthogonality  of E : Mean absolute error %8.3e',...
		 checksit1))
    difer(checksit1,[],[],NaN)
  end
  
  % Compare the sum over the eigenvalues with the partial Shannon number
  Nm=indeks(nsubm(N,m,1),m+1);
  disp(sprintf('Nm = %4.2g ; V1 = %4.2e',Nm,V(1)))

  % Only check if you've calculated at least 2N of them
  if K>=max(2*N,3)
    difer(sum(V)-Nm,8,1,'Partial Shannon numbers check out')
  end
    
  if ~isnan(x)
    % Make them start on an upswing if at all possible - might be avoided
    % if we wanted to keep the eigenvalue signed, which makes not much
    % sense of course, given the quadratic concentration problem
    for index=1:size(E,2)
      E(:,index)=E(:,index)*sign(indeks(E(~~E(:,index),index),1));
      E2(:,index)=E2(:,index)*sign(indeks(E2(~~E2(:,index),index),1));
    end
    
    % And remember, if you didn't like the scaling, you can undo it here
    if scalem==0 && ~strcmp(method,'GL')
      warning off MATLAB:divideByZero
      E=E./repmat(sqrt(x(:)),1,size(E,2));
      E2=E2./repmat(sqrt(x(:)),1,size(E2,2));
      disp('First element will be NaN and the FFT will be bogus')
      warning on MATLAB:divideByZero
    end
  end
  
  % Keep track of computation time
  ts=etime(clock,t0);
  disp(sprintf('\nComputation took %8.4f s',ts))

  % Warn if things are getting tight
  if any(V<eps)
    warning('Numerical degeneracy suspected of near-zero eigenvalues')
  end
  
  
  
  % Provide output
  vars={E,V,Nm,c,x,E2};
  varargout=vars(1:nargout); 
elseif strcmp(m,'demo1')
  % DeVilliers et al. (2003) Figure 1 using two methods
  ah(1)=subplot('position',[0.13 0.4096 0.775 0.5154]);
  m=0;
  x=linspace(0,1.5,200);
  [E,V,Nm,c]=swdisk(m,25,4,[],x,'DV');
  plot(x,E(:,1:4),'Marker','o','LineS','none','MarkerS',4);
  l=legend('1','2','3','4','Location','SouthWest');
  hold on
  [E1,V,Nm,c]=swdisk(m,25,4,[],x,'SE');
  plot(x,E1(:,1:4),'LineW',2); 
  axis tight; ylim([-1 1]); plot([1 1],ylim,'k-')
  xl(1)=xlabel('x'); ylabel('\phi(x)'); title(sprintf('c = %g ; m = %g',c,m))
  grid on; hold off
  ah(2)=subplot(3,1,3);
  plot(x,E-E1); hold on
  axis tight; ylim([-1 1]*1e-14); plot([1 1],ylim,'k-')
  grid on; hold off
  xl(2)=xlabel('x'); ylabel('\phi_{DV}(x)-\phi_{SE}(x)')
  title(sprintf('Misfit between different methods'))
  longticks(ah,2); delete(xl(1)); nolabels(ah(1),1)
  fig2print(gcf,'portrait')
  figdisp
  % Now list the eigenvalues
  format long ; sprintf('%.7e\n',V)
  format short
elseif strcmp(m,'demo2')
  % Slepian (1964), Figure 3a
  clf
  m=0; N=1; method='GL';
  [E,V,Nm,c,x]=swdisk(m,N,4,[],linspace(0,1.5,100),method);
  % Do the normalization differently from our own - to 1 over the interval
  E=sqrt(2*pi)*E(:,1:4)./sqrt(repmat(V,size(E,1),1));
  plot(x,E,'LineW',2);
  l=legend('1','2','3','4','Location','SouthWest');
  axis([0 1.5 -4 4])
  xlabel('x') 
  ylabel('\phi(x)')
  title(sprintf('m = %i ; c = %g ; N = %g ; meth = %s',m,c,N,method))
  set(gca,'xtick',[0:0.1:1.5]); hold on
  plot([1 1],[-4 7],'k')
  grid on; hold off
  fig2print(gcf,'portrait'); 
  figdisp([],2)
elseif strcmp(m,'demo3')
  % Slepian (1964), Figure 3b
  clf
  m=0; N=25; method='SE';
  [E,V,Nm,c,x]=swdisk(m,N,4,[],linspace(0,1.5,100),method);
  % Do the normalization differently from our own - to 1 over the interval
  E=sqrt(2*pi)*E(:,1:4)./sqrt(repmat(V,size(E,1),1));
  plot(x,E,'LineW',2); 
  l=legend('1','2','3','4','Location','SouthWest');
  axis([0 1.5 -4 4])
  xlabel('x') 
  ylabel('\phi(x)')
  title(sprintf('m = %i ; c = %g ; N = %g ; meth = %s',m,c,N,method))
  set(gca,'xtick',[0:0.1:1.5]); hold on
  plot([1 1],[-4 7],'k')
  grid on; hold off
  fig2print(gcf,'portrait'); 
  figdisp([],3)
elseif strcmp(m,'demo4')
  % Slepian (1964), Figure 4a
  clf
  m=2; N=1; method='GL';
  [E,V,Nm,c,x]=swdisk(m,N,4,[],linspace(0,1.5,100),method);
  % Do the normalization differently from our own - to 1 over the interval
  E=sqrt(2*pi)*E(:,1:4)./sqrt(repmat(V,size(E,1),1));
  plot(x,E,'LineW',2); 
  l=legend('1','2','3','4','Location','SouthWest');
  axis([0 1.5 -4 4])
  xlabel('x') 
  ylabel('\phi(x)')
  title(sprintf('m = %i ; c = %g ; N = %g ; meth = %s',m,c,N,method))
  set(gca,'xtick',[0:0.1:1.5]); hold on
  plot([1 1],[-4 7],'k')
  grid on; hold off
  fig2print(gcf,'portrait'); 
  figdisp([],5)
elseif strcmp(m,'demo5')
  % Slepian (1964), Figure 4b
  clf
  m=2; N=25; method='GL';
  [E,V,Nm,c,x]=swdisk(m,N,4,[],linspace(0,1.5,100),method);
  % Do the normalization differently from our own - to 1 over the interval
  E=sqrt(2*pi)*E(:,1:4)./sqrt(repmat(V,size(E,1),1));
  plot(x,E,'LineW',2); 
  l=legend('1','2','3','4');
  axis([0 1.5 -4 7])
  xlabel('x') 
  ylabel('\phi(x)')
  title(sprintf('m = %i ; c = %g ; N = %g ; meth = %s',m,c,N,method))
  set(gca,'xtick',[0:0.1:1.5]); hold on
  plot([1 1],[-4 7],'k')
  grid on; hold off
  fig2print(gcf,'portrait'); 
  figdisp([],4)
elseif strcmp(m,'demo6')
  % Shkolnisky (2007), Figure 1
  N=(1/2).^2;
  m=1;
  ens=[0 1 5 10];
  methods={'DV','GL','SE'};
  method=methods{ceil(abs(3*guess(1))/10)};
  method='DV';
  [E,V,Nm,c,x,E2]=swdisk(m,N,max(ens)+1,[],[],method);
  clf
  plot(x,E2(:,ens+1),'LineW',2);
  for index=1:length(ens)
    legs{index}=sprintf('%s = %i','\alpha',ens(index)+1);
  end
  l=legend(legs,'Location','SouthWest');
  axis([0 1 -6 8])
  xlabel('x')
  ylabel('\phi(x)')
  title(sprintf('m = %i ; c = %g ; N = %g ; meth = %s',m,c,N,method))
  grid on
  fig2print(gcf,'portrait')
elseif strcmp(m,'demo7')
  % Shkolnisky (2007), Figure 2
  N=(10/2).^2;
  ems=[0 1 5 10]; 
  methods={'DV','GL','SE'};
  method=methods{ceil(abs(3*guess(1))/10)};
  [a,V,Nm,c,x]=swdisk(ems(1),N,2,[],[],method);
  E(:,1)=a(:,2)/sqrt(V(2))*sqrt(2*pi);
  for index=2:length(ems)
    [a,V]=swdisk(ems(index),N,2,[],[],method);
    E(:,index)=a(:,2)/sqrt(V(2))*sqrt(2*pi);
  end
  for index=1:length(ems)
    legs{index}=sprintf('%s = %i','m',ems(index));
  end
  clf
  plot(x,E(:,1:length(ems)),'LineW',2);
  l=legend(legs,'Location','SouthWest');
  axis([0 1 -5 3])
  xlabel('x')
  ylabel('\phi(x)')
  title(sprintf('Method %s',method))
  grid on
  fig2print(gcf,'portrait')
elseif strcmp(m,'demo8')
  % Shkolnisky (2007), Figure 3
  N=([0.1 1 10 100]/2).^2;
  m=1; 
  methods={'DV','SE'};
  method=methods{ceil(abs(2*guess(1))/10)};
  [a,V,Nm,c,x,b]=swdisk(m,N(1),2,[],[],method);
  E(:,1)=b(:,2);
  for index=2:length(N)
    [a,V,Nm,c,x,b]=swdisk(m,N(index),2,[],[],method);
    E(:,index)=b(:,2);
  end
  clf
  plot(x,E(:,1:length(N)),'LineW',2);
  l=legend('c = 0.1','c = 1','c = 10','c = 100');
  axis([0 1 -3 3])
  xlabel('x')
  ylabel('\phi(x)')
  title(sprintf('Method %s',method))
  grid on
  fig2print(gcf,'portrait')
elseif strcmp(m,'demo9')
  % Slepian (1964), Table 1
  N=([0.1 0.5 1 1.5 2 3 4 5 10]/2).^2;

  m1=0;
  for index=1:length(N)
    [E,V1(index,:)]=swdisk(m1,N(index),1,[],'SE');
  end
  
  N=([0.5 1 2 3 4 5 6 10]/2).^2;
  m2=1;
  for index=1:length(N)
    [E,V2(index,:)]=swdisk(m2,N(index),1,[],'SE');
  end

  N=([1:10]/2).^2;
  m3=0;
  for index=1:length(N)
    [E,V3(index,:)]=swdisk(m3,N(index),2,[],'SE');
  end

  N=([1 2 4:10]/2).^2;
  m4=1;
  for index=1:length(N)
    [E,V4(index,:)]=swdisk(m4,N(index),2,[],'SE');
  end

  N=([1 2 5:16]/2).^2;
  m5=0;
  for index=1:length(N)
    [E,V5(index,:)]=swdisk(m5,N(index),3,[],'SE');
  end

  N=([1 2 5:14]/2).^2;
  m6=1;
  for index=1:length(N)
    [E,V6(index,:)]=swdisk(m6,N(index),3,[],'SE');
  end

  N=([1 2 5 7:16]/2).^2;
  m7=0;
  for index=1:length(N)
    [E,V7(index,:)]=swdisk(m7,N(index),4,[],'SE');
  end

  N=([1 2 5 9:17]/2).^2;
  m8=1;
  for index=1:length(N)
    [E,V8(index,:)]=swdisk(m8,N(index),4,[],'SE');
  end

  N=([1:10]/2).^2;
  m9=2;
  for index=1:length(N)
    [E,V9(index,:)]=swdisk(m9,N(index),1,[],'SE');
  end

  N=([1 2 5:11]/2).^2;
  m10=2;
  for index=1:length(N)
    [E,V10(index,:)]=swdisk(m10,N(index),2,[],'SE');
  end

  N=([1 2 5 7:15]/2).^2;
  m11=2;
  for index=1:length(N)
    [E,V11(index,:)]=swdisk(m11,N(index),3,[],'SE');
  end

  N=([1 2 5 10:18]/2).^2;
  m12=2;
  for index=1:length(N)
    [E,V12(index,:)]=swdisk(m12,N(index),4,[],'SE');
  end

  disp(sprintf('\nm = %i ; alpha = %i\n',m1,1))
  disp(sprintf('%.7e\n',V1))

  disp(sprintf('\nm = %i ; alpha = %i\n',m2,1))
  disp(sprintf('%.7e\n',V2))

  disp(sprintf('\nm = %i ; alpha = %i\n',m3,2))
  disp(sprintf('%.7e\n',V3(:,2)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m4,2))
  disp(sprintf('%.7e\n',V4(:,2)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m5,3))
  disp(sprintf('%.7e\n',V5(:,3)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m6,3))
  disp(sprintf('%.7e\n',V6(:,3)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m7,4))
  disp(sprintf('%.7e\n',V7(:,4)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m8,4))
  disp(sprintf('%.7e\n',V8(:,4)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m9,1))
  disp(sprintf('%.7e\n',V9(:,1)))
  
  disp(sprintf('\nm = %i ; alpha = %i\n',m10,2))
  disp(sprintf('%.7e\n',V10(:,2)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m11,3))
  disp(sprintf('%.7e\n',V11(:,3)))

  disp(sprintf('\nm = %i ; alpha = %i\n',m12,4))
  disp(sprintf('%.7e\n',V12(:,4)))
elseif strcmp(m,'demo10')
  % Compare the three methods
  m=abs(guess(1));
  N=abs(guess(1)*2)+10;
  K=abs(guess(1)); K=K+(K==0);
    
  format long
  
  x=linspace(0,1.5,100);
  [E1,V1,Nm1,c1,x1]=swdisk(m,N,K,[],x,'DV');
  [E2,V2,Nm2,c2,x2]=swdisk(m,N,K,[],x,'SE');
  [E3,V3,Nm3,c3,x3]=swdisk(m,N,K,[],x,'GL');

  % Check eigenvalues
  difer(V1-V2)
  difer(V2-V3)
  difer(V1-V3)
  % Check eigenfunctions where they're supposed to be good
  difer(E1(x<=1,:)-E2(x<=1,:))
  difer(E2(x<=1,:)-E3(x<=1,:))
  difer(E1(x<=1,:)-E3(x<=1,:))
  clf
  a=plot(x,E1,'k+'); hold on
  b=plot(x,E2,'-'); 
  c=plot(x,E3,'o'); 
  hold off
  ylim(1.1*minmax(E2(x<=1,:)))
  set(gca,'xtick',[0 1 2],'xgrid','on')
  %aha=nanmean(E3./E2);
  %plot(1:K,aha,'o-')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,P2]=getscaledjacobis(x,m,L)
x=x(:)';
P=nan(L+1,length(x));
P2=P;
% Do the next thing by recursion all-at-once later
% Inkling was that the normalization can be absorbed into the DeVilliers'
% matrix such that it becomes symmetric as we expect from a self-adjoint
% operator. Discussion with Yoel Shkolnisky on 07/11/2012.
for l=0:L
  P(l+1,:)=x.^(m+1/2).*jacobi(l,m,1-2*x.^2)*...
	   factorial(m)/prod(l+[1:m]);
  P2(l+1,:)=x.^(m+1/2).*jacobi(l,m,1-2*x.^2)*sqrt(2*(2*l+m+1));
end

% Check this makes sense
if ~~sum(isnan(P(:)))
  error('Trouble - decrease L or switch method to ''SE''')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,J2]=getscaledbessels(x,m,L,c)
x=x(:)';
J=[besselj(m+2*[0:L]+1,c*x(:))]';
J2=J;
warning off MATLAB:divideByZero
for l=0:L
  J(l+1,:) =J(l+1,:)./prod(m+[1:l])*factorial(l)./sqrt(c*x);
  J2(l+1,:)=J2(l+1,:).*sqrt(2*(2*l+m+1))./sqrt(c*x);
end
warning on MATLAB:divideByZero
% Put in the limit form from A.S. 9.1.7 and l'Hopital's rule
J(:,x==0)=0;
J2(:,x==0)=0;


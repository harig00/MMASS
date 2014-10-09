function varargout=svdslep(N,NW,K,method,imp)
% [E,V,SE,ngro]=SVDSLEP(N,NW,K,method,imp)
%
% Explicit diagonalization of the Slepian concentration operators 
% in one linear dimension.  
%
% INPUT:
%
% N       The length of the sequence [default: 2^8]
% NW      Half the Shannon number (time-bandwidth product) [default: 10] 
% K       The requested number of eigenvectors [default: 2NW-1]
% method  1 Using EIGS on the Hermitian form [default]
%         2 Using EIG on the Hermitian form
%         3 Using SVDS on the projection operator
%         4 Using SVD on the projection operator
%         5 Using the power method for the largest eigenvalue
% imp     1 Uses the implicit operator method [default]
%         0 Uses the sparse explicit operator method [slow]
%
% OUTPUT:
%
% E     The eigenfunctions of the concentration problem
% V     The eigenvalues of the concentration problem
% SE    The power spectrum of the eigenfunctions
% ngro  The "growth" factor used
%
% EXAMPLE:
%
% svdslep('demo1')
% 
% Last modified by fjsimons-at-alum.mit.edu, 06/21/2010

% See some notes in RB VIII p 90

% Default values
defval('N',2^10)

if ~isstr(N)
  defval('NW',10)
  defval('K',2*NW-1)
  defval('method',1)
  defval('ngro',256);
  defval('imp',1);

  % This one for "excessive" verification
  defval('xver',0)

  % Make a larger domain
  Nd=N*ngro;

  if ngro>32 || N>256
    % Force to the implicit method or your machine will die
    imp=1;
    disp('Method reset to IMPLICIT')
  end
  
  if imp==1
    disp(sprintf('Implicit method embedded %ix\n',ngro))
    nonz=(Nd-N)/2+[1:N];
    P =@(x) proj(x,nonz);
    Q =@(x) fft(x);
    Qi=@(x) ifft(x);
    nonz=[1:ngro*NW+1 Nd-ngro*NW+1:Nd];
    L =@(y) proj(y,nonz);
    H =@(x) P(Qi(L(Q(P(x)))));
    % Acknowledge that H is complex (though it is symmetric)
    OPTS.isreal=false;
    OPTS.disp=0;
    % Remember to specify the output size
    [E,V]=eigs(H,Nd,K,'LR',OPTS);
    [V,i]=sort(diag(V),'descend');
    E=E(:,i); V=V(1:K); E=E(:,1:K);
  elseif imp==0
    disp(sprintf('Explicit, sparse method embedded %ix\n',ngro))
    % The spatial projection operator
    P=diag([zeros(1,(Nd-N)/2) ones(1,N) zeros(1,(Nd-N)/2)]);
    % The Fourier transform operator
    Q=dftmtx(Nd)/sqrt(Nd);
    % The bandlimiting operator
    L=diag([ones(1,ngro*NW+1) zeros(1,Nd-2*ngro*NW-1) ones(1,ngro*NW)]);

    % Make them all sparse for speed - could start from it altogether
    P=sparse(P);
    Q=sparse(Q);
    L=sparse(L);
    
    if xver==1
      % Check that Q is a unitary matrix
      difer(Q'*Q-diag(diag(Q'*Q)))
      difer(Q'*Q-eye(size(Q)))
      % Check that H is a Hermitian matrix
      difer(H'-H)
      % Check that H is a positive definite matrix
      [R,p]=chol(H); isposdef=p==0; difer(isposdef-1)
      % Check that Q does what I think it does
      f=rand(Nd,1);
      % See the normalization of FFT (1) and IFFT (1/N)
      difer(fft(f)/sqrt(Nd)-Q*f);
    end
    
    % The operator whose singular functions we want
    A=L*Q*P;
    
    % The operator whose eigenfunctions we want
    H=P*Q'*A;
    
    % The singular values of A are the eigenvalues of AA'. The singular
    % values of the concentration/projection operator are the eigenvalues of
    % the "squared" projection operator.

    % The eigenvector decomposition
    switch method
     case 1
      disp('Using EIGS')
      OPTS.disp=0;
      [E,V]=eigs(H,K,'LR',OPTS);
      [V,i]=sort(diag(V),'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 2   
      disp('Using EIG')
      [E,V]=eig(H,'nobalance');
      [V,i]=sort(diag(V),'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 3
      disp('Using SVDS')
      [U,V,E]=svds(A,K);
      [V,i]=sort(diag(V).^2,'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 4
      disp('Using SVD')
      [U,V,E]=svd(A);
      [V,i]=sort(diag(V).^2,'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 5
      disp('Using the power method for the first one only')
      Vi=1; Vj=2;
      E=rand(ngro*N,1);
      E=E/sqrt(E'*E);
      while abs(Vj-Vi)>1e-6
	Vi=Vj;
	Vj=E;
	E=H*E;
	Vj=real(E'*Vj);
	E=real(E/sqrt(E'*E));
	plot(E); 
	pause(0.1)
      end
      V=Vj;
    end
  end  
  
  % Define some kind of tolerance level
  tol=100*eps;

  % Make them real as we know they should be
  if any(imag(V)>tol)
    error('Complex eigenvalues');
  else
    V=real(V);
    % Check imaginary part of the "good" eigenfunctions
    disp(sprintf('mean(abs(imag(E))) = %8.3e out of %8.3e\n',...
		 mean(mean(abs(imag(E(:,V>tol))))),...
		 mean(mean(abs(E(:,V>tol))))))

    % Take out only the central part
    E=E((Nd-N)/2+1:(Nd-N)/2+N,:);

    % Note that they were normalized in the complex plane
    E=real(E); E=E./repmat(diag(sqrt(E'*E))',size(E,1),1);
    % Cannot do abs(E).*sign(E) because the sign is bad for complex
    % But either real or imag seems to work... they more or less scale!
  end

  if xver==1
    % Check the orthogonality for the "good" eigenfunctions
    ortho=E(:,V>tol)'*E(:,V>tol);
    difer(ortho-diag(diag(ortho)))
    difer(diag(ortho)-1)
  end

  if nargout>2
    % Get the power spectrum
    SE=abs(fft(E,8*length(E))).^2;
  else
    SE=NaN;
  end
  
  % Fix the first bit to be an upswing
  for index=1:size(E,2)
    E(:,index)=E(:,index)*sign(E(2,index));
  end
  
  % Output
  varns={E,V,SE,ngro};
  varargout=varns(1:nargout);
elseif strcmp(N,'demo1')
  % So you can say svdslep('demo1',1)
  defval('NW',[])
  imp=NW;
  defval('imp',1)
    
  N=2^(3+ceil(rand*5)); NW=ceil(rand*5);
  N=2^8; NW=3;
  method=1;
  disp(sprintf('\nN = %i , NW = %i\n',N,NW))
  [E,V,SE,ngro]=svdslep(N,NW,max(3,2*NW),method,imp);
  % Compare with Matlab's built-in functions
  [E2,V2]=dpss(N,NW,max(3,2*NW));
  
  clf
  [ah,ha]=krijetem(subnum(2,3));

  % Make the frequency axis for the positive frequencies only
  [fax,selekt]=fftaxis1d(E(:,1),length(SE),N-1);
  SE=SE(selekt,:);

  % Another, identical, more recently revised, way
  % [K,kx]=knum2([2 N],[2 N]-1);
  % fx=-fliplr(indeks(kx,selekt)/2/pi);

  yls=[-100 2];

  for ind=1:3
    SE2=indeks(abs(fft(E2(:,ind),8*length(E2)).^2),selekt);

    axes(ah(ind))
    plot(E(:,ind),'b','linew',2); hold on
    plot(E2(:,ind),'y')
    
    disp(sprintf('Eigenfunction %i: rms %8.3f%s',ind,...
		 rms(E(:,ind)-E2(:,ind))/rms(E(:,ind))*100,'%'))
    
    yl(ind)=ylabel('Slepian eigenfunction');

    set(ah(ind),'xtick',[1 N],'xlim',[1 N])
    
    xlabel(sprintf('%s = %12.9f (EIGS)','\lambda',V(ind)))
    
    axes(ha(2*ind))
    plot(N*fax,decibel(SE(:,ind)),'b','linew',2); hold on
    plot(N*fax,decibel(SE2),'y');
    plot(N*[NW/N NW/N],yls,'k-');
    plot(-N*[NW/N NW/N],yls,'k-');
    hold off

    %xl(3+ind)=xlabel('frequency (sample)^{-1}');
    xl(3+ind)=xlabel('frequency x sample size');
    yl(3+ind)=ylabel('power spectral density (dB)');
    title(sprintf('%s = %12.9f (DPSS)','\lambda',V2(ind)))

    % Show the entire range
    too=8*1.95;
    exes=[0 1/too];
    exes=[-0.1 20];
    set(ha(2*ind),'xtick',round(exes*10)/10,'xlim',minmax(exes),...
		  'xgrid','off','ylim',yls)
  end
  
  % Cosmetics
  delete(yl([2 3 5 6]))
  longticks(ah)
  
  tt=supertit(ah(1:3),...
	   sprintf('SVDSLEP vs DPSS on a %ix domain',ngro),12);
  movev(tt,0.15)
  fig2print(gcf,'landscape')
  
  % Is the ratio of the two things symmetric?
  %   for ind=1:size(E,2)
  %     a=E(:,ind)./E2(:,ind);
  %     rmsym(ind)=rms([a(1:length(a)/2)-flipud(a(length(a)/2+1:end))])./...
  % 	  rms(a)*100;
  %   end
  %   rmsym

  figdisp([],'demo1')  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pv=proj(v,p)
% Projects the vector v on the indices p
Pv=zeros(size(v));
Pv(p)=v(p);

function varargout=svdslep2(N,R,K,method,imp)
% [E,V,SE]=SVDSLEP2(N,R,K,method,imp)
%
% Explicit diagonalization of the Slepian concentration operator
% in two dimensions.  
%
% INPUT:
%
% N       The length of one side of the domain [default: 2^4]
% R       The Shannon ratio (of full frequency spectrum) [default: .1]
% K       The number of eigenvectors [default: 10]
% method  1 Using EIGS on the Hermitian form [default]
%         2 Using EIG on the Hermitian form
%         3 Using SVDS on the projection operator
%         4 Using SVD on the projection operator
%         5 Using the power method for the largest eigenvalue
% imp     1 Uses the implicit operator method [default]
%         0 Doesn't 
%
% OUTPUT:
%
% E     The eigenfunctions of the concentration problem
% V     The eigenvalues of the concentration problem
% SE    The power spectrum of the eigenfunctions
%
% EXAMPLE:
%
% svdslep2('demo1')
%
% SEE ALSO: SVDSLEP3
%
% REMARKS: Stability remains a problem... even the results between the
% sparse and the non-sparse explicit approaches are different in the
% details (though with identical eigenvalues and with vast computation
% speed differences.)
%
% Written by Eugene Brevdo and Frederik J Simons, 04/03/2010
% Last modified by fjsimons-at-alum.mit.edu, 05/27/2011

% Default values
defval('N',2^6)

if ~isstr(N)
  defval('R',.1);
  defval('K', 10);
  defval('method',1)
  defval('ngro',8);
  defval('imp',1);

  % This one for "excessive" verification
  defval('xver',0)

  % Make a larger domain, inflate one size
  Nd=N*ngro;

  disp(sprintf('ngro = %i\n',ngro))

  if ngro>sqrt(32) || N>sqrt(256)
    % Force to the implicit method or your machine will die
    imp=1;
    disp('Method reset to IMPLICIT')
  end
  
  % Make the two-dimensional spatial projection operator...
  nonz=matranges(...
      repmat((Nd-N)/2*(Nd+1)+[1 N],N,1)'+repmat([0:Nd:Nd*(N-1)]',1,2)');
  if imp==0
    [nonzi,nonzj]=ind2sub([Nd^2 Nd^2],nonz+(nonz-1)*Nd^2);
    P=sparse(nonzi,nonzj,1,Nd^2,Nd^2);
    disp(sprintf('Explicit, sparse method embedded %ix\n',ngro))
    % Projection operator built from the image mask
    if xver==1
      PI=zeros(Nd); 
      % ... see CENTERVAL
      PI((Nd-N)/2+1:(Nd-N)/2+N,(Nd-N)/2+1:(Nd-N)/2+N)=1;
      difer(nonz(:)-find(PI))
      difer(P-diag(PI(:)));
      PI=proj(ones(Nd),nonz);
      difer(P-diag(PI(:))); 
    end
  elseif imp==1
    P=@(x) proj(x,nonz);
    disp(sprintf('Implicit method embedded %ix\n',ngro))
  end
  
  % Make the two-dimensional spectral projection operator...
  if imp==0
    % The Fourier transform operator
    Q=sparse(dft2mtx(Nd)/Nd);
  elseif imp==1
    Q= @(x) fft2vec(x);
    Qi=@(y) ifft2vec(y);
  end

  % The bandlimiting operator
  LI=zeros(Nd);

  % Calculate distance from origin (0,0)
  [LIx,LIy]=meshgrid(linspace(-1/sqrt(2),1/sqrt(2),Nd));
  LIx=fftshift(LIx); LIy=fftshift(LIy);
  % Find the elements within radius R of (0,0) (all are within R=1)
  LI=sqrt(LIx.^2 + LIy.^2)<=R;
  % Draw a box of side length percentage R
  % LI=(sqrt(2)*abs(LIx)<=R) & (sqrt(2)*abs(LIy)<=R);

  if imp==0
    % Projection operator built from image mask
    L=sparse(diag(LI(:)));
    
    % The operator whose singular functions we want
    A=L*Q*P;
    
    % The operator whose eigenfunctions we want
    H=P*Q'*A;
    % The singular values of A are the eigenvalues of AA'. The singular
    % values of the concentration/projection operator are the eigenvalues of
    % the "squared" projection operator.

    if xver==1
      % Check that Q is a unitary matrix
      difer(Q'*Q-diag(diag(Q'*Q)))
      difer(Q'*Q-eye(size(Q)))

      % Check that H is a Hermitian matrix
      difer(H'-H)
      % Check that H is a positive definite matrix
      [R,p]=chol(H); isposdef=p==0; difer(isposdef-1)
      % Check that Q does what I think it does
      f=rand(Nd,Nd);
      % See the normalization of FFT (1) and IFFT (1/N)
      difer(fft2(f)/Nd-reshape(Q*f(:),Nd,Nd));
    end
  else
    % But how does it know that the output is 2D?
    L=@(y) proj(y,find(LI));
    H=@(x) P(Qi(L(Q(P(x)))));
    % Test if you want: H(rand(Nd*Nd,1))
  end

  if imp==1
    % Acknowledge that H is complex (though it is symmetric)
    OPTS.isreal=false;
    OPTS.disp=0;
    % Remember to specify the output size
    [E,V]=eigs(H,Nd^2,K,'LR',OPTS);
    
    [V,i]=sort(diag(V),'descend');
    E=E(:,i); V=V(1:K); E=E(:,1:K);
  elseif imp==0
    % The eigenvector decomposition
    switch method
     case 1
      OPTS.disp=0;
      [E,V]=eigs(H,K,'LR',OPTS);
      [V,i]=sort(diag(V),'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 2   
      [E,V]=eig(H,'nobalance');
      [V,i]=sort(diag(V),'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 3
      [U,V,E]=svds(A,K);
      [V,i]=sort(diag(V).^2,'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 4
      [U,V,E]=svd(A);
      [V,i]=sort(diag(V).^2,'descend');
      E=E(:,i); V=V(1:K); E=E(:,1:K);
     case 5
      Vi=1; Vj=2;
      E=rand(N,1);
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
  tol=sqrt(eps); %100*eps;

  % Make them real as we know they should be
  if any(imag(V)>tol)
    error('Complex eigenvalues');
  else
    V=real(V);
    % Check imaginary part of the "good" eigenfunctions
    disp(sprintf('mean(abs(imag(E))) = %8.3e out of %8.3e\n',...
		 mean(mean(abs(imag(E(:,V>tol))))),...
		 mean(mean(abs(E(:,V>tol))))))
    % Note that they were normalized in the complex plane
    E=real(E); E=E./repmat(diag(sqrt(E'*E))',size(E,1),1);
    % Should this be a signed ABS rather than a real?
  end
  
  % Take out only the part we care about: the central image, size N
  % Note that the solution is actually LIMITED to the central part
  E=E(nonz, :);

  if xver==1
    % Check the orthogonality for the "good" eigenfunctions
    ortho=E(:,V>tol)'*E(:,V>tol);
    difer(diag(ortho)-1)
    difer(ortho-diag(diag(ortho)))
    % If it doesn't, sometimes, it could be due to numerical degeneracy on
    % the eigenvalues. Let's just live with it.
  end

  if nargout>2
    % Get the power spectrum
    SE=zeros(N^2,size(E,2));
    for i=1:size(E,2),
      SE(:,i)=indeks((fftshift(abs(fft2(v2s(E(:,i)))).^2)),':');
    end
  else
    SE=NaN;
  end

  % Output
  varns={E,V,SE,ngro};
  varargout=varns(1:nargout);
elseif strcmp(N,'demo1')
  % So you can say svdslep('demo1',1)
  defval('R',[])
  imp=R;
  defval('imp',1)
  
  N=2^6;
  [E,V,SE,ngro]=svdslep2(N,[],[],[],imp);

  % Check these out - with imp=0 get some funny degeneracies even for the
  % square case!
  V
  
  clf
  [ah,ha]=krijetem(subnum(2,3));

  for ind=1:3
    axes(ah(ind))
    imagesc(v2s(E(:,ind)))
    axis image
    set(ah(ind),'xtick',[1 N],'xlim',[1 N],...
		'ytick',[1 N],'ylim',[1 N])
    xlabel(sprintf('%s = %12.9f (EIGS)','\lambda',V(ind)))

    % Plot the power spectral density
    axes(ha(2*ind))
    imagesc(decibel(v2s(SE(:,ind))))
    axis image

  end
  
  % Cosmetics
  longticks(ah)
  fig2print(gcf,'landscape')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pv=proj(v,p)
% Projects the vector v on the indices p
Pv=zeros(size(v));
Pv(p)=v(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fv=fft2vec(v)
% Returns the two-dimensional FFT of a vector
Fv=fft2(v2s(v));
Fv=Fv(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iFv=ifft2vec(Fv)
% Returns the two-dimensional IFFT of a vector
iFv=ifft2(v2s(Fv));
iFv=iFv(:);


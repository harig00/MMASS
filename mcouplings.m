function M=mcouplings(L,Lmax,approx,appL)
% M=mcouplings(L,Lmax,approx,appL)
%
% Calculates the mode-coupling matrix for eigenvalue-weighted multitaper
% analysis of arbitrary regions. Its shape depends only on bandwidth L.
% Dahlen and Simons (2008) eq. (145).
%
% INPUT 
%
% L        Bandwidth (maximum angular degree) of the taper
% Lmax     Bandwidth (maximum angular degree) of the spectrum
% approx   0 No approximations are made, the whole thing is calculated
%          1 When the second index exceeds appL, the rows are copied down
% appL     From where on is the approximation invoked [default: 3*L]
%          Note: will need to calculate (appL+L) Wigner 0j coefficients
%
% OUTPUT:
%
% M        Coupling matrix (the eigenvalue-weighted sum of the
%          individual-taper coupling matrices MG from MCOUPLING)
%
% See also: MCOUPLING
%
% Last modified by fjsimons-at-alum.mit.edu, 04/02/2007

% Maybe later modify to do only certain rows

defval('L',18)
defval('xver',0)
defval('Lmax',100)
defval('approx',0)

if approx==1
  defval('appL',3*L)
  if appL+L>=Lmax
    % It's not an approximation now anymore is it
    approx=0;
    appL=0;
    disp('You''re getting an exact result')
  end
else
  appL=0;
end

if approx~=0 & approx~=1
  error('Specify valid option approx')
end

fname=fullfile(getenv('IFILES'),'MCOUPLINGS',...
	       sprintf('MCOUPLINGS-%i-%i-%i-%i.mat',...
		       L,Lmax,approx,appL));

if exist(fname)==2
  disp(sprintf('load %s',fname))
  load(fname)
else
  tic
  % Initialize matrix with zeros
  M=repmat(0,[Lmax+1 Lmax+1]);
  
  h=waitbar(0,'MCOUPLINGS: Doing the sums, be patient');

  LmaxT=Lmax;
  if approx==1
    Lmax=appL;
  end
  
  % Best to get all of the Wigner 0-j symbols at once
  % This is always faster than computing unsaved 0j's on the fly
  [jk,C0,S0]=zeroj(0,0,0,Lmax+L*(approx==1));
  
  % Do the lower triangular half of the coupling matrix
  for l=0:Lmax
    % Know this is a BANDLIMITED kernel so just stop at the right time
    % This cuts the computation time dramatically
    for lp=max(0,l-L):l+L*(approx==1)
      % Calculate the arrays of Wigner 0j symbols
      %FJS WAG=zeroj([0:L],repmat(l,1,L+1),repmat(lp,1,L+1),...
      WAG=zeroj([0:L],l,lp,...
		Lmax+L*(approx==1),[],C0,S0).^2*(2*[0:L]+1)';
      if xver==1
	% This is the old way where you compute these on the fly
	WAG=wigner0j(L,l,lp).^2*(2*[0:L]+1)';
	difer(WAG-WAG2,[],[],'MCOUPLINGS Check passed')
	%disp('Another excessive test passed')
      end
      % Fill the elements
      M(l+1,lp+1)=WAG;
    end
    waitbar(lp/(Lmax+1),h)
  end
  
  if approx==0
    % Symmetrize this portion appropriately in one line
    M=[M-diag(diag(M))]+[M-diag(diag(M))]'+diag(diag(M));
    % How about: M=M+tril(M',-1);
  end
  
  % Don't forget about the (2l'+1) for the column dimensions
  % M=M.*repmat(2*[0:LmaxT]+1,LmaxT+1,1)/(L+1)^2;
  % This is the same thing, perhaps a bit more elegantly
  % Remember it's constant diagonal AFTER the asymmetry
  M=M*diag(2*[0:LmaxT]+1)/(L+1)^2;

  if approx==1
    % Use an approximation
    % The row index of the last row
    i=Lmax+1;
    % The column index of the last row
    j=Lmax-L+1;
    % Extract the nonzero elements from the last row and copy it down
    s=repmat(M(i,j:j+2*L),LmaxT-Lmax,1);
    % Make column indices that repeat    
    j=pauli(j+1:LmaxT+1+L,2*L+1);
    i=repmat([i+1:LmaxT+1]',1,2*L+1);
    % This is the continuation of the lower triangular part
    S=sparse(i,j,s);
    M=M+S(:,1:LmaxT+1);
  end

  % Save for future reference
  eval(sprintf('save %s M L Lmax LmaxT',fname))
  toc
  close(h)
end

% Check we're dealing with the right normalization etc.
difer(sum(M(1:LmaxT-L+1,:),2)-1,9,[],...
      sprintf('MCOUPLINGS: Row sum equals %5.3f check passed',...
	      mean(sum(M(1:LmaxT-L+1,:),2))))



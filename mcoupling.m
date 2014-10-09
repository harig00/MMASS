function [MG,V,N,MH]=mcoupling(TH,L,Lmax,sord)
% [MG,V,N,MH]=mcoupling(TH,L,Lmax,sord)
%
% Calculates the mode-coupling matrix for single- and multi-taper
% analysis of axisymmetric single or double polar caps. The estimated
% spectrum out to degree Lmax is the true spectrum linearly transformed
% by this (Lmax+1)X(Lmax+1) matrix, which contains the effect of a
% degree-L taper. Dahlen and Simons 2008 eq. (131).
%
% INPUT 
%
% TH       Angular extent of the spherical cap, in degrees OR
% L        Bandwidth (maximum angular degree) of the taper
% Lmax     Bandwidth (maximum angular degree) of the spectrum
% sord     1 Single polar cap of diameter TH [default]
%          2 Double polar cap of diameter 2TH
%
% OUTPUT:
%
% MG       A three-dimensional array with the (L+1)^2 single-taper
%          couplings; the number of tapers is the third dimension; for
%          the bandlimited tapers
% V        The sorted eigenvalues belonging to every one of the tapers
% N        The Shannon number
% MH       Same, but for the space-limited tapers
%
% EXAMPLE:
%
% MG=mcoupling; tn=1;
% for l=0:100 ; plot([0:100],MG(l+1,:,tn),'-o'); 
% grid on; ylim(minmax(MG(:,:,tn))); longticks(gca,2)
% set(gca,'xtick',unique([0 max(0,l-18) l min(l+18,100) 100])); 
% pause; end
%
% See also: MCOUPLINGS
%
% Last modified by fjsimons-at-alum.mit.edu, 05/16/2007

defval('TH',10)
defval('L',18)
defval('Lmax',100)
defval('xver',1)
defval('sord',1)

fname=fullfile(getenv('IFILES'),'MCOUPLING',...
	       sprintf('MCOUPLING-%i-%i-%i-%i.mat',TH,L,Lmax,sord));

if exist(fname)==2
  disp(sprintf('load %s',fname))
  load(fname)
else
  % First get the sum over all orders of the square of the coefficients
  % of all tapers, collected per taper (A) and per order (L) in GAL
  [G,V,EL,EM,N,GAL]=glmalpha(TH,L,sord);
  % You'll want to sort this
  [V,i]=sort(V,'descend');
  GAL=GAL(i,:);
  if nargout>3
    [H,V,EL,EM,N,HAL]=hlmalpha(TH,L,2*Lmax,sord);
    HAL=HAL(i,:);
    MH=repmat(NaN,[Lmax+1 Lmax+1 (L+1)^2]);
  end
    
  % Initialize matrix
  MG=repmat(0,[Lmax+1 Lmax+1 (L+1)^2]);
  
  h=waitbar(0,'MCOUPLING: Doing the sums, be patient');

  % Do the upper triangular half of the coupling matrix
  for l=0:Lmax
    for lp=max(0,l-L):l
      % Don't divide by 4pi since we've already done that in GLMALPHA by
      % normalizing to one; instead of by normalizing to 4pi and then
      % dividing. The column dimension is the taper number now 
      % Note that the GAL already contains the (2l+1), i.e. we really
      % need the total power and not the power spectral density
      WAG=wigner0j(L,l,lp).^2*GAL';
      % To check the unbiasedness...
      % Note that (2*[0:L]+1).*wigner0j(L,l,lp) only sums to one if
      % L=l+lp; but sum(2*[0:2*Lmax]+1).*wigner0j(2*Lmax,l,lp)==1

      % Fill the elements; this part is symmetric (the kernel is NOT)
      MG(l+1,lp+1,:)=WAG;
      MG(lp+1,l+1,:)=WAG;

      if nargout>3
	WAH=wigner0j(2*Lmax,l,lp).^2*HAL';
	MH(l+1,lp+1,:)=WAH;
	MH(lp+1,l+1,:)=WAH;
      end
    end
    waitbar(lp/(Lmax+1),h)
  end
  close(h)
  % Don't forget about the (2l'+1) for the column dimensions
  MG=MG.*repmat(2*[0:Lmax]+1,[Lmax+1 1 (L+1)^2]);
  if nargout>3
    MH=MH.*repmat(2*[0:Lmax]+1,[Lmax+1 1 (L+1)^2]);
  else
    MH=NaN;
  end
  eval(sprintf('save %s MG MH TH L Lmax V N',fname))
end

% Excessive verification
if xver==1
  M1=zeros(Lmax+1);
  for index=1:length(V)
    M1=M1+MG(:,:,index)*V(index)/sum(V);
  end
  % Compare with the analytic formula of the eigenvalue-weighted sum
  difer(M1-mcouplings(L,Lmax),[],[],'MCOUPLING: Check passed');
  % Should perhaps verify that the sum(M,2) is one until Lmax-L?
end

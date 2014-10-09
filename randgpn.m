function varargout=randgpn(k,dc,dcn,xver)
% [Zk,Zx]=randgpn(k,dc,dcn,xver)
% 
% Returns a set of (complex proper) normal variables suitable for
% IFFT, most notably this works for even and/or odd sized rectangles.
%
% INPUT:
%
% k       A wavenumber matrix (DC component in center)
% dc      The (m,n) indices to the DC component in k
% dcn     The (m,n) indices to the Nyquist components in k
%         --> Note that these three are straight out of KNUM2, and that
%         the actual wavenumbers are not used, only size(k) is needed. 
% xver    1 Checks the Hermiticity of the result by inverse transformation
%
% OUTPUT:
%
% Zk      A size(k) matrix with (complex proper) normal Fourier variables
% Zx      A size(k) matrix with a real-valued random field
%
% EXAMPLE 1
%
% [k,kx,ky,dci,dcn]=knum2([round(rand*100) round(rand*100)],[100 100]);
% [Zk,Zx]=randgpn(k,dci,dcn); imagesc(Zx)
%
% EXAMPLE 2 
% 
% [k,kx,ky,dci,dcn]=knum2([200 201],[100 100]);
% F=abs(randn*10); [Zk,Zx]=randgpn(k,dci,dcn); Zf=ifft2(ifftshift(Zk.*exp(-F*k)));
% imagef([0 0],size(Zf),Zf); axis image ;
% title(sprintf('exp(%3.3gk)',-F)); colorbar('hor')
%
% SEE ALSO: KNUM2
%
% Last modified by fjsimons-at-alum.mit.edu, 08/12/2013

defval('xver',0)

% Make a receptacle with the dimension of the wavenumbers
[m,n]=size(k);
Zk=zeros(m,n); clear k

% Determine the parity
modd=mod(m,2);
nodd=mod(n,2);

% Define the ranges to pick out a half plane; DC is (0,0) is included
lhn=1:dc(2);
rhn=dc(2)+1:n;

% Now make some complex proper normal random variables of half variance
ReZk=randn(m,dc(2))/sqrt(2);
ImZk=randn(m,dc(2))/sqrt(2);
% And make a handful real normal random variables of unit variance
realZk=randn(4,1);

% Fill the entire left half plane with these random numbers
lh=ReZk+sqrt(-1)*ImZk;
Zk(:,lhn)=lh;

% Now overwrite and symmetrize the column containing zero
Zk(dc(1)+1:end,dc(2))=conj(flipud(Zk(2-modd:dc(1)-1,dc(2))));

% Fill the center with a real, where the DC component goes 
Zk(dc(1),dc(2))=realZk(1);

% Fill the Nyquist with a real, if we capture them exactly 
for ind=1:size(dcn,1)
  Zk(dcn(ind,1),dcn(ind,2))=realZk(1+ind);
end

% The intersection of two columns with a Nyquist should also be real
% cos(k_x x+k_y y)+i sin(k_x x+k_y y) at (k_x,k_y)=(-pi,-pi) and for x and
% y integers, the sin term vanishes thus the cos is real and the
% multiplier must be real! And the x and y are integers because they are
% multiplying the sampling interval which factors out. 
Zk(~modd,~nodd)=realZk(1+ind+1);

% And now symmetrize the right half plane so that output is Hermitian
Zk(2-modd:end,rhn)=conj(flipud(fliplr(lh(2-modd:end,2-nodd:dc(2)-1))));

% Symmetrize the top row that may contain the Nyquist, if it does
Zk(~modd,dc(2)+1:end)=conj(fliplr(Zk(~modd,2-nodd:dc(2)-1)));

% Symmetrize the leftmost column that may contain the Nyquist, if it does
Zk(dc(1)+1:end,~nodd)=conj(flipud(Zk(2-modd:dc(1)-1,~nodd)));

if nargout>1 || xver==1
  Zx=ifft2(ifftshift(Zk));
  % You could now check that the IFFT2 is real (don't forget the "2"!!)
  if ~isreal(Zx) ; error(sprintf('Not Hermitian by %5g',mean(imag(Zx(:))))); end
else
  Zx=NaN;
end

% Output
varns={Zk,Zx};
varargout=varns(1:nargout);

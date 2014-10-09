function [xout,a,d,err,PA,PS]=fbank(x,tipe,N,M,pph)
% [xout,a,d,err,PA,PS]=fbank(x,tipe,N,M,pph)
%
% Combines ABANK and SBANK so the perfect reconstruction
% condition becomes immediately testable.
%
% 'pph' 1= Time-domain full bitrate (inefficient);
%       2= Time-domain polyphase (inefficient);  
%       3= Z-domain polyphase (default)
%
% Returns the Frequency-domain (pph=3) 
%          or Time domain (pph=3) polyphase matrix as
%          PA (P-analysis) and PS (P-synthesis)
%
% Single-iteration filter-bank of signal 'x' with filters
% of a certain 'name'.
%
% Returns reconstructed signal without delay or distortion,
% and the rms error of the approximation.
%
% See also WTLIFT.

% Last modified by fjsimons-at-alum.mit.edu, 03/01/2002

t0=clock;

defval('pph',3)

[h0,f0]=wc(tipe,[N M]);

% Check normalization of the product filter and
% create synthesis coefficients from the analysis ones:
% conjugate mirror filters with alternating signs to satisfy
% the no-alias condition
[h0,h1,f0,f1,l,lf0,lh0]=prodco(f0,h0);
% Get downsampling arrays based on the data length
% and the filter lengths.
[d2fx,d2hx,lx,lxout]=landd(f0,h0,x);

% Get polyphase matrices as cell arrays and in components
[PA,PS,I,...
 H0even,H0odd,H1even,H1odd,...
 F0even,F0odd,F1even,F1odd]=polyphase(f0,h0);

% Analysis
switch pph
 case 3
  %  This is the "Standard Algorithm"
  disp('Z-domain polyphase implementation')
  Xeven=x(even(x));
  Xodd=x(~even(x));
  % This here is the polyphase matrix in the z-form (not Toeplitz):
  % Our definition is the transpose of Daubechies'
  % It's a block matrix with Laurent polynomial coefficients  
  % Always need to offset them, but, depending on the length of the  
  % input series, may need to add an extra zero or not.
  % If h0 and f0 are both odd, get different numbers of evens
  % and odds.  
  % We're taking care of the case were lf0 and lh0 are both
  % even and thus split into equal lengh parts
  lh0lxe=length(H0even)+length(Xeven);
  lh0lxo=length(H0odd)+length(Xodd);
  lh1lxe=length(H1even)+length(Xeven);
  lh1lxo=length(H1odd)+length(Xodd);
  
  % This is tricky business, but it works
  % -> for even- and odd length signals
  % -> for even- and odd length filters
  a=repmat(0,max(lh0lxe,lh0lxo-1),1);
  d=repmat(0,max(lh1lxe,lh1lxo-1),1);
  
  a(2:lh0lxe)=conv(H0even(:),Xeven(:));  
  a(1:lh0lxo-1)=a(1:lh0lxo-1)+conv(H0odd(:),Xodd(:));

  d(2:lh1lxe)=conv(H1even(:),Xeven(:));
  d(1:lh1lxo-1)=d(1:lh1lxo-1)+conv(H1odd(:),Xodd(:));

  case 2
    disp('Time-domain polyphase implementation')
    % Filter Toeplitz matrices
    % Size of H0 is (lx+lh0-1) by (lx)
    % Size of H1 is (lx+lf0-1) by (lx)
    % and element (1,1) is h0(1)
    H0=convmtx(h0(:),lx);
    H1=convmtx(h1(:),lx);
    % Downsample: delete rows: keep these
    AL=H0(d2hx,:);
    AB=H1(d2fx,:);
    % Polyphase form
    % Split signal in even and odd
    % Note that this is dependent on your first sample:
    % do you call it x(0) or x(1)
    % The Matlab convention differs from the others,
    % so our even is their odd.
    X=[x(even(x)) ; x(~even(x))];
    % Construct the polyphase ANALYSIS matrix (Toeplitz)
    H0even=AL(:,even(x));
    H0odd=AL(:,~even(x));
    H1even=AB(:,even(x));
    H1odd=AB(:,~even(x));
    % Implement filters by multiplication
    % This is the polyphase applied to the samples in the time domain
    PA=[H0even H0odd ;
        H1even H1odd];
    ad=PA*X;
    % This is the intermediate vector; split up if you want
    a=ad(1:sum(d2hx));
    d=ad(sum(d2hx)+1:end);
  case 1
    disp('Toeplitz implementation')
    % Filter Toeplitz matrices
    % Size of H0 is (lx+lh0-1) by (lx)
    % Size of H1 is (lx+lf0-1) by (lx)
    % and element (1,1) is h0(1)
    H0=convmtx(h0(:),lx);
    H1=convmtx(h1(:),lx);
    % Downsample: delete rows: keep these
    AL=H0(d2hx,:);
    AB=H1(d2fx,:);  
    % Apply - this yields the decomposition we're after
    a=AL*x(:);
    d=AB*x(:);
    % And they have length ceil((lx+lh0-1)/2) and ceil((lx+lf0-1)/2)
    % Polyphase matrix is not applicable
    PA=[];
 otherwise
  error('Specify valid method')
end

switch pph
 case 3
  lf0ea=length(F0even)+length(a);
  lf1ed=length(F1even)+length(d);
  lf0oa=length(F0odd)+length(a);
  lf1od=length(F1odd)+length(d);  

  % This is tricky business, but it works
  % -> for even- and odd length signals
  % -> for even- and odd length filters
  y0=repmat(0,floor(lxout/2),1);
  y1=repmat(0,ceil(lxout/2),1);

  xout=repmat(NaN,lxout,1);
  
  y0(1:lf0ea-1)=conv(F0even(:),a(:));
  y0(1:lf1ed-1)=y0(1:lf1ed-1)+conv(F1even(:),d(:));

  y1(1:lf0oa-1)=conv(F0odd(:),a(:));
  y1(1:lf1od-1)=y1(1:lf1od-1)+conv(F1odd(:),d(:));
     
  % Even output
  xout(even(xout),1)=y0;
  % Odd output
  xout(~even(xout),1)=y1;
    
  case 2
  % Synthesis
  % Size of F0 is (lx+lh0-1+lf0-1) by (lx+lh0-1)
  % Size of F1 is (lx+lf0-1+lh0-1) by (lx+lf0-1)
  % and element (1,1) is h0(1)
  % Filter Toeplitz matrices
  F0=convmtx(f0(:),lx+lh0-1);
  F1=convmtx(f1(:),lx+lf0-1);
  % Construct the polyphase SYNTHESIS matrix
  % Polyphase form
  % Upsample: delete columns: keep these
  F0=F0(:,d2hx);
  F1=F1(:,d2fx);
  % This already reduces the number of multiplications by two
  % Now split in odd and even
  xout=repmat(NaN,lxout,1);
  F0even=F0(even(xout),:);
  F0odd=F0(~even(xout),:);
  F1even=F1(even(xout),:);
  F1odd=F1(~even(xout),:);
  PS=[F0even F1even ;
      F0odd F1odd];
  % Could make this TYPE II (p. 132)
  % The output of this is mixed and it needs to be recombined!  
  % Check out the operation PS*PA*X - you've switched evens 
  % and odds and added a delay here and there so the 
  % reconstruction automatically makes up for it
  xoutm=PS*ad;
  % Output of the original evens of x
  xout(even(xout),1)=xoutm(1:size(F0even,1));
  % Output of the original odds of x
  xout(~even(xout),1)=xoutm(size(F0even,1)+1:end);
case 1
  % Synthesis
  % Size of F0 is (lx+lh0-1+lf0-1) by (lx+lh0-1)
  % Size of F1 is (lx+lf0-1+lh0-1) by (lx+lf0-1)
  % and element (1,1) is h0(1)
  % Filter Toeplitz matrices
  F0=convmtx(f0(:),lx+lh0-1);
  F1=convmtx(f1(:),lx+lf0-1);
  % Upsample
  blk0=~(1:lx+lh0-1);
  blk1=~(1:lx+lf0-1);
  blk0(d2hx)=a;
  blk1(d2fx)=d;
  w1=F0*blk0(:);  
  w2=F1*blk1(:);
  xout=w1+w2;
  % Polyphase is not applicable
  PS=[];
 otherwise
  error('Specify valid method')
end

% Complete filter bank; length is shifted by the power of
% the only odd power of z^-1 that has a nonzero coefficient,
% i.e. by (length(P0)-1)/2
% Can also do wkeep(xout,size(x)).
xout=xout(l+1:end-l);

err=sqrt(mean((xout(:)-x(:)).^2));

if err>1e-10
  error(sprintf('Reconstruction failed with %.2f',err))
  else
    disp(sprintf('FBANK Reconstruction error %8.4f',err))
end
disp(sprintf('FBANK (Analysis & Synthesis) took %8.4f s',etime(clock,t0)))

% Given the condition that we will choose 'conjugate mirror filters'
% F_0(z)=+H_1(-z)
% F_1(z)=-H_0(-z)
% to automatically satisfy the 'no-aliasing' requirement (4.2),
% the 'no-distortion' requirement (4.3) leads to finding
% P_0(z)=F_0(z)H_0(z)=z^-l where l is odd and even powers of z
% can be multiplied by any coefficient. We choose a lowpass
% P_0(z)=(1+z^-1)^2pQ(z) where Q(p) is a polynomial of degree 2p-2.
% Next, we can factor P_0 into different choices of H_0 and F_0,
% for various time lags l, and retrieve F_1 and H_1 from there.
% We take symmetric product filters, so l leads to 2l roots.
% The choice for l=1 is the Haar bank, l=3 leads to a 6th degree
% polynomial which can be split into two third degree polynomials,
% hence two 4-tap filters, called Daubechies-4. These are causal,
% but not linear phase.

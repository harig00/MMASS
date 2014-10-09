function H=hermitian(mat)
% H=hermitian(mat)
%
% Transforms a complex vector of even length or a square complex matrix
% with sides of even length into a Hermitian one.
%
% ODD LENGTH, says SCO, is best. EVEN, if improperly done at the Nyquist,
% will lead to funky behavior that we see.
% 
% INPUT: 
%
% A Fourier matrix after FFTSHIFT (long wavelengths in center)
%
% OUTPUT:
%
% A Hermitian matrix, after IFFTSHIFT, so you can proceed straight to
% IFFT, which will be real.
%
% Hermitian symmetry is when
% w[n]=w*[-n], where w* denotes the complex conjugate of w.
% If x element of Real^N, then X is Hermitian.
% Care is taken of special points.
%
% See KNUM2, SYNCOH, SYNSPEC, ACTSPEC and CPYSPEC
%
% Example:
%
% [fld,sspec,k,ftherm]=synspec(128);
% subplot(211); plot(real(ftherm))
% subplot(212); plot(imag(ftherm))
% for index=1:length(ftherm)/2
%    verify(index)=...
%          ftherm(index+1)==...
%          conj(ftherm(length(ftherm)-index+1));
% end
% So the FIRST sample is the odd one out, and 
% length(ftherm)/2+1 is compared with itself.
%
% The action of FFTSHIFT is to put the longest wavelength in the center.
%
% ftherm=fftshift(ftherm);
% for index=length(ftherm)/2+2:length(ftherm)
%    verify(index-length(ftherm)/2-1)=...
%             ftherm(index)==...
%             conj(ftherm(length(ftherm)-index+2));
% end
% so  length(ftherm)/2+1 is the odd one out,
% and 1 is compared with itself.
%
% Note that when plotting this using fftaxis, need to FLIPLR!
% Because in FFTAXIS, floor((dim+1)/2) is the zero-frequency sample.
% In a way, it would be easier to have an odd sequence length.
%
% Frederik Simons with Oded Aharonson
% Last modified by fjsimons-at-mit.edu, May 20th, 2002

if any(size(mat)==[1 1])  
  % For vectors
  if any(mod(length(mat),2))
    error('Input length must be mod-2')
  end  
  n=length(mat);
  nh=n/2;
  nh1=nh+1;
  mat(1)=real(mat(1));
  mat(nh1:n)=conj(flipud(fliplr(mat(2:nh1))));
  mat(~isfinite(mat))=0; 
  H=ifftshift(mat);
elseif size(mat,1)==size(mat,2)
  % For matrices
  if any(mod(size(mat),2))
    error('Input length must be mod-2')
  end
  n=size(mat,1);
  nh=n/2;
  nh1=nh+1;
  % Lower right quadrant = conj transpose of upper left quadrant
  mat(nh1:n,nh1:n)=conj(flipud(fliplr(mat(2:nh1,2:nh1))));
  % Lower left quadrant = conj transpose of upper right quadrant
  mat(nh1:n,2:nh1)=conj(flipud(fliplr(mat(2:nh1,nh1:n))));
  % The first rowof UR = first row of UL
  mat(1,nh1:n)=conj(flipud(fliplr(mat(1,2:nh1))));
  % First colum of LL = first column of UL
  mat(nh1:n,1)=conj(flipud(fliplr(mat(2:nh1,1))));
  % The DC term is invariant
  mat(nh1,nh1)=real(mat(nh1,nh1));
  % The special point in the middle is invariant
  mat(1,1)=real(mat(1,1));
  % The halfway point on the first row is invariant
  mat(1,nh1)=real(mat(1,nh1));
  % The halfway point in the first column is invariatint
  mat(nh1,1)=real(mat(nh1,1));
  % Everything must be finite or else set to zero
  mat(~isfinite(mat))=0;
  % The first element must be the minimum of the whole thing
  % And then you can shift things over
  mat(1,1)=min(real(mat(:)));
  H=ifftshift(mat);  
else
  error('Not a time series or a square matrix')
end






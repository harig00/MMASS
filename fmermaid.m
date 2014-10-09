function x=fmermaid(x)
% x=FMERMAID(x)
%
% Forward MERMAID-transform.
%
% INPUT:     
%
% x             The integer real-valued signal of length a power of two
%
% OUTPUT:
%
% x             The same length Mallat-ordered integer wavelet transform
%
% COMMENTS:
%
% Performs 5 iterations of a CDF(2,4) lifting wavelet on the input
% signal which consists of a power of two number of samples. All
% hardwired. No boundary filters are used in this case, and no scaling
% constants are included. Yet. See ABANK, SBANK, and MAKEWC on how to do
% this properly, which is by translating the intermediate scalings into
% lifting steps also.
%
% The transform is linear: difer(2*fmermaid(x)-fmermaid(x*2))
% 
% See also IMERMAID, FMERMAIDTEST
%
% Last modified by fjsimons-at-alum.mit.edu, 18.04.2006

% Random input values
defval('x',round(randn(1024,1)*round(rand*18)));

% Check the values are integers to begin with
if ~all(round(x)==x)
  error('Input series must consist of integer values')
end

% Check the length
if 2^nextpow2(length(x))~=length(x)
  error('Input array must have length a power of two')
end

% Input is a column vector
x=x(:);

% Initial length
lx=length(x);

% Number of cascades
ncasc=5;

% Loop over the filter bank branches as a cascading filter bank
for index=1:ncasc
  % PREDICT step to make the lowpass/approximation/scaling coefficients
  % Replace the even samples with their prediction from the odd ones
  for n=2:2:lx-2
    x(n)=x(n)-floor([x(n-1)+x(n+1)]/2+1/2); 
  end
  % UPDATE step to make the highpass/detail/wavelet coefficients
  % Update the odd samples based on the new even samples
  for n=5:2:lx-3
    x(n)=x(n)+floor([-3*x(n-3)+19*x(n-1)+19*x(n+1)-3*x(n+3)]/64+1/2);
  end
  % REARRANGE to a Mallat multiresolution organization, not in-place
  % First the approximation, then the wavelet coefficients
  x=[x(1:2:lx) ; x(2:2:lx) ; x(lx+1:end)];
  % REFINE: 
  % In the next iteration, work with the approximation coefficients
  lx=length(x)/2^index;
end











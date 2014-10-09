function x=imermaid(x)
% x=IMERMAID(x)
%
% Inverse MERMAID-transform.
%
% INPUT:     
%
% x       The Mallat-ordered integer wavelet transform of length a power
%         of two 
%
% OUTPUT:
%
% x       The same length integer real-valued signal 
%
% COMMENTS:
%
% Simply forget about the constant here, let the terms grow,
% you have an invertible transform... add them in later if you 
% want, or take them into account when you need to compare them.
% If not, while it may sometimes work, it actually usually doesn't... 
%
% See also FMERMAID.
%
% Last modified by fjsimons-at-alum.mit.edu, 18.04.2006

% Check the length
if 2^nextpow2(length(x))~=length(x)
  error('Input array must have length a power of two ')
end

% Work with a column vector 
x=x(:);

% Number of cascades
ncasc=5;

% Loop over the filter bank branches as a cascading filter bank
for index=ncasc:-1:1
  % Increasing length considered
  lx=length(x)/2^index;
  % UNDO UPDATE step to make new odds
  for n=3:lx-1
    x(n)=x(n)-floor([-3*x(lx+n-2)+19*x(lx+n-1)+...
		     19*x(lx+n)-3*x(lx+n+1)]/64+1/2);
  end
  % UNDO PREDICT step to make new evens
  for n=1:lx-1
    x(lx+n)=x(lx+n)+floor([x(n)+x(n+1)]/2+1/2); 
  end  
  % REARRANGE
  % These are the values we just made
  y0=x(1:lx);
  % And here is the next batch of the wavelet coefficients
  y1=x(lx+1:lx*2);
  % The odd values get to be here
  x(1:2:lx*2)=y0;
  % The even values get to be here
  x(2:2:lx*2)=y1;
end

































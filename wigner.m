function varargout=wigner(sig,sint,tax)
% WIGNER(sig,sint)
% WIGNER(sig,sint,tax)
% [abstf,fax,tax]=WIGNER(...)
%
% 'sig' is the signal
% 'sint' is the sampling interval
% 'tax' is a time axis [DEFAULT linspace(0,nsig*sint,nsig)]
% 'fax' is the frequency axis
%
% Computes the alias-free Wigner-Ville Distribution
%
% See also WVD, WVDist, WVDist_AF

% Ref: Jeong and Williams, IEEE-SP 40 (11) 2757-2765

% Written by fjsimons-at-mit.edu, October 3rd, 2000

sig = sig(:);
if length(sig)~=2.^nextpow2(length(sig)); error('Dyadic length required');  end

nsig   = length(sig);
nfft=256;

f   = [zeros(nsig,1); sig; zeros(nsig,1)];
afwig = zeros(nfft, nsig);
ix  = 0:(nsig/2-1);
zerosn = zeros(nsig,1);

for t=1:nsig,
  tplus    = nsig + t + ix;
  tminus   = nsig + t - ix;
  x = zerosn;
  % Even indices
  x(1:2:nsig) = f(tplus) .* f(tminus);
  % Odd indices
  x(2:2:nsig) = (f(tplus+1).*f(tminus) + f(tplus).*f(tminus-1))/2;  
  afwig(:, t) = 2* (fft(x,nfft));
end

abstf = real(afwig);

% Make frequency and time axes
[fax,selekt]=fftaxis1D(sig,nfft,(nsig-1)*sint);
if nargin==2
  tax=linspace(0,(nsig-1)*sint,nsig);
end
abstselekt=abstf(selekt,:);

if nargout
  nam={ 'abstselekt','fax','tax'};
  for index=1:length(nam);
    varargout{index}=eval(nam{index});
  end
else
  contourf(tax,2*pi*fax,abstselekt);
  colorbar('ver')
  axis('xy')
  title('Alias Free Wigner Distribution');
  xlabel('Time')
  ylabel('Angular Frequency \omega')
end


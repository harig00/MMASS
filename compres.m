function [cf,rmse,acor,ener]=compres(x,xs,an,dn)
% [cf,rmse,acor,ener]=COMPRES(x,xs,an,dn)
%
% Figures out compression and quality
% of wavelet transforms.
% Works with the result of 
% SUCCAPP or IWT
%
% INPUT:
%
% x       Original signal
% xs      Matrix with successive approximations
% an      Number of scaling coefficients
% ds      Numbers of wavelet coefficients
%
% OUTPUT:
% 
% cf      Compression factor (per cent, cumulatively per scale)
% rmse    Relative RMS Error of the approximation (0-100%)
% acor    Autocorrelation peak of the approximation (0-100%)
% ener    In terms of the energy of the approximation (0-100%)
%
% Last modified by fjsimons-at-alum.mit.edu, 07/14/2008

nuco=cumsum([an(end) fliplr(dn)]);

cf=100*nuco/length(x);
cxs=cumsum(xs,2);

if nargout>1
  for index=1:size(xs,2)
    rmse(index)=100-(rms(x(:)-cxs(:,index))/rms(x(:))*100);
    acor(index)=100*max(xcorr(x(:),cxs(:,index),'coeff'));
    ener(index)=100*rms(cxs(:,index))^2/rms(x(:))^2;
  end
end








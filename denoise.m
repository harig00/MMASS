function varargout=denoise(x,tipe,nvm,n)
% [xr,cf,rmse,acor,ener]=denoise(x,tipe,nvm,n)
% 
% Denoises a signal by wavelet thresholding
%
% INPUT
%
% x       Input signal
% tipe    'Daubechies' or 'CDF'
% nvm     Number of vanishing moments 
% n       Number of scales used in decomposition
%
% OUTPUT
%
% xr      Reconstructed signal after thresholding
% cf      Compression factor
% rmse    RMSE error
% acor    Cross-correlation error
% ener    Energy error
%
% EXAMPLES:
%
% denoise('demo1')
%
% SEE ALSO:
%
% THRESHOLD, THRESHOLD2
%
% Last modified by fjsimons-at-alum.mit.edu, 05/27/2009

defval('tipe','Daubechies')
defval('nvm',4)
defval('n',nextpow2(length(x))-1)

if ~isstr(x)
  % Calculate by polyphase since there we have the 
  % most extensive library
  pph=3;
  % Perform wavelet transform
  [a,d,an,dn,ts,cf]=wt(x,tipe,nvm,n,pph);
  % Perform thresholded inverse wavelet transform - see also THRESHOLD2
  [dt,dnz]=threshold(d,dn,'soft');
  [xc,xr]=iwt(a,dt,an,dn,tipe,nvm);
  if nargout>1
    cf=100*sum([an(end) dnz])/length(x);
    rmse=100-(rms(x(:)-xr(:))/rms(x(:))*100);
    acor=100*max(xcorr(x(:),xr(:),'coeff'));
    ener=100*rms(xr(:))^2/rms(x(:))^2;
  end
  varns={'xr' 'cf' 'rmse' 'acor' 'ener'};
  for index=1:nargout
    varargout{index}=eval(varns{index});
  end
elseif strmatch(x,'demo1')
  t=0.001:0.001:2; 
  y=chirp(t,0,1,10);       
  yn=sigmerge(y,randn(size(y)));
  subplot(211)
  plot(t,yn)
  hold on
  plot(t,y,'y','linew',2)
  hold off
  subplot(212)
  xr=denoise(yn,'Daubechies',4,5);
  plot(t,yn)
  hold on  
  plot(t,sum([xr{:}],2),'y','linew',2)
end

function varargout=syncoh(field1,coh2)
% field2=SYNCOH(field1,coh2)
%
% Calculates a second field with a known coherence with respect to a
% first field.
%
% INPUT:
%
% field1       A matrix with some kind of real-valued data
% coh2         A matrix with some kind of coherence function
% 
% OUTPUT:
%
% field2       A matrix with the second field of real-valued data
%
% EXAMPLE: 
%
% syncoh('demo1')
% syncoh('demo2')
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2008

defval('field1','demo1')

if ~isstr(field1)
  [irow,icol]=size(field1);
  F=fftshift(fft2(field1-mean(field1(:))));
  var=1;
  
  % Construct real and imaginary part of noise
  rand1=var*randphase(irow,icol);
  rand2=var*randphase(irow,icol);
  
  % Protect coherence from blowing up
  coh2(coh2<0.01)=0.01;
  
  % Lowry and Smith, 1994, Eq. 7
  absn2=abs(F).^2.*(1-coh2)./coh2;
  
  ra2=abs(rand1).*absn2;
  ra=sqrt(ra2).*sign(rand1);
  ia=abs(sqrt(absn2-ra2).*sign(rand2));
  
  rg=ra+real(F);
  ig=ia+imag(F);
  
  % Make Hermitian - this is key
  H=hermitian(rg+i*ig);
  field2=real(ifft2(H));
  
  % Produce output
  varns={field2};
  varargout=varns(1:nargout);
elseif strcmp(field1,'demo1')
  ddir= '/home/fjsimons/MyPapers/2003/JGR-2003/DATA/COHSYNTHETIC';
  load(fullfile(ddir,'Xmidcpy'))
  field1=topocpy;
  
  load(fullfile(ddir,'cohdata'))
  
  field2=syncoh(field1,C2K);

  subplot(221)
  imagesc(field1); axis image
  subplot(222)
  imagesc(C2K); axis image
  subplot(223)
  imagesc(field2); axis image; 
  caxis(prctile(field2(:),[2.5 97.5]))
  subplot(224)
  [FX,FY,SX,SY,SXY,COH2]=mtm3(field1,field2,[],[]);
  imagesc(COH2); axis image
elseif strcmp(field1,'demo2')
  x=linspace(-10,10,108); y=x;
  [XR,YR]=rotdom(x,y,pi/6);
  coh2=scale(gauss2(XR,YR,2,4),[0 1]);

  ddir= '/home/fjsimons/MyPapers/2003/JGR-2003/DATA/COHSYNTHETIC';
  load(fullfile(ddir,'Xmidcpy'))
  field1=topocpy;
  
  load(fullfile(ddir,'cohdata'))
  
  field2=syncoh(field1,coh2);
  
  subplot(221)
  imagesc(field1); axis image
  subplot(222)
  imagesc(coh2); axis image
  subplot(223)
  imagesc(field2); axis image; 
  caxis(prctile(field2(:),[2.5 97.5]))
  subplot(224)
  [FX,FY,SX,SY,SXY,COH2]=mtm3(field1,field2,[],[]);
  imagesc(COH2); axis image
end



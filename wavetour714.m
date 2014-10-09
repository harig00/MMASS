function wavetour714
% WAVETOUR714
%
% Reproduces Figure 7.14 of Mallat's "Wavelet Tour" book.
% using my own codes.
%
% Last modified by fjsimons-at-alum.mit.edu, 21.04.2006

N=[2 4 ; 3 7];

[ah,ha]=krijetem(subnum(2,2*size(N,1)));

for index=1:size(N,1)
  [phix,phi,wx,w,philim,wlim]=graphs('CDF',N(index,1:2));

  [phixd,phid,wxd,wd,philimd,wlimd]=graphs('CDF',N(index,1:2),1);

  axes(ha(4*index-3))
  pp=plot(phix,phi,'LineW',1);
  axis tight

  axes(ha(4*index-1))
  pp=plot(phixd,phid,'LineW',1);
  axis tight

  axes(ha(4*index-2))
  pw=plot(wx,w,'LineW',1);
  axis tight
  
  axes(ha(4*index))
  pw=plot(wxd,wd,'LineW',1);
  axis tight
end

% Cosmetics
set(ah,'ylim',[-1.5 2])
fig2print(gcf,'portrait')
figdisp

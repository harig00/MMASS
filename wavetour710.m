function wavetour710
% WAVETOUR710
%
% Reproduces Figure 7.10 of Mallat's "Wavelet Tour" book.
% using my own codes.
%
% Last modified by fjsimons-at-alum.mit.edu, 21.04.2006

N=[2 3 4];

[ah,ha]=krijetem(subnum(2,length(N)));

for index=1:length(N)
  [phix,phi,wx,w,philim,wlim]=graphs('Daubechies',N(index));
  axes(ah(index))
  pp=plot(phix,phi,'LineW',1);
  axis tight

  axes(ah(index+length(N)))
  pw=plot(wx,w,'LineW',1);

  axis tight
end

% Cosmetics
set(ah,'ylim',[-1.5 2])
fig2print(gcf,'portrait')
figdisp

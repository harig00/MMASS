function tenlectures63
% TENLECTURES63
%
% Reproduces Figure 6.3 of Daubechies "Ten Lectures" book.
% using my own codes.
%
% Last modified by fjsimons-at-alum.mit.edu, 21.04.2006

N=[2 3:2:9];

[ah,ha]=krijetem(subnum(length(N),2));

for index=1:length(N)
  [phix,phi,wx,w,philim,wlim]=graphs('Daubechies',N(index));
  axes(ha(index))
  pp=plot(phix,phi,'LineW',1);
  axis tight

  axes(ha(index+length(N)))
  pw=plot(wx,w,'LineW',1);

  axis tight
end

% Cosmetics
set(ah,'ylim',[-1.5 2])
fig2print(gcf,'tall')
figdisp

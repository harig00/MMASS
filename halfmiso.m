function [khalf,g2M,knumi,ravi]=halfmiso(data,kxc,kyc)

% [khalf,g2M,knumi,ravi]=halfmiso(data,kxc,kyc)
%
% For a two-dimensional matrix 'data' calculates the location of the maximum
% value along the isotropic centered average as well as the location of the
% half maximum. For symmetric square matrices. 
% Put in the complete (neg and pos) wavenumber axis 'kxc' and 'kyc'.
% We assume this is a symmetric matrix with identical wave numbered ranges in
% the x and y direction.
%
% Return half wavenumber 'khalf' (for g2M/2) and maximum coherence 'g2M'.
% 
% See also HALFMAZI, the azimuthal equivalent.

% Last modified by fjsimons-at-mit.edu, November 19th, 2001

precirat=1/40;

[ry,rx]=size(data);
hrx=ceil(rx/2);
hry=ceil(ry/2);
%if hrx~=hry
%  error('Should input square data set')
%end

[knum,rav]=isav(data,kxc,kyc);

if sum(isnan(rav))>0
  error('NaNs in average will mess up interpolation') 
end

[aj,callpro]=star69;

knumi=linspace(knum(1),knum(end),10*length(knum));
ravi=interp1(knum,rav,knumi,'cubic');
[g2M,kMi]=max(ravi);
% Not the zero-wavelength
kMi=kMi+1*(kMi==1);

preci=g2M*precirat;
howclose=abs(ravi-g2M/2);
k12Mi=min(find(howclose<preci));
% NOT NEXT LINE; WANT FIRST WAVENUMBER THAT COMPLIES!
% [mi,k12Mi]=min(howclose);

if strcmp(callpro,'australiadir3') | ...
      strcmp(callpro,'australiadir40') 
  hold on
%  plot(knum,rav,'r+');
  plot(knumi,ravi,'g-','LineW',2); hold on ; set(gca,'xscale','log')
  xlim([min(kxc(kxc>0)) kxc(end)])
  plot(knumi(kMi),g2M,'o','MarkerF','g','MarkerE','k')
  plot(knumi(k12Mi),g2M/2,'o','MarkerF','g','MarkerE','k')
  xlo=xlim;
%  plot(xlo,[g2M/2 g2M/2],'--');
%  plot(xlo,[g2M g2M],'--','Color','k'); 
%  plot(xlo,[g2M/2 g2M/2]+preci,'--','Color',[1 1 1]*0.7); 
%  plot(xlo,[g2M/2 g2M/2]-preci,'--','Color',[1 1 1]*0.7); 
  ylim([0 1])
  longticks
end

khalf=knumi(k12Mi);

if isempty(khalf); khalf=0; end


function varargout=cohen(funname,varargin)
% COHEN('hermite',lwin,R,[])
% COHEN('hermite',lwin,R,nwin)
% COHEN('slepian',lwin,TW,TC,nwin)
% [totabs,fax,tax]=COHEN(...)
%
% Calculates the smoothing kernel for the Cohen's class
% interpretation of a spectrum estimation
% Calculates the average Wigner-Ville spectrum
% for a set of functions, weighted by their
% eigenvalues, which depend on the radius R of the 
% circular time-frequency domain.
%
% See also HERMITE, LAGUERRE, BAYRAM
%
% EXAMPLE I
%
% cohen('hermite',128,3,6)
%
% EXAMPLE II (Lilly and Park, Fig. 2)
%
% lwin=128; nwin=6; TW=3; TC=8.0;
% cohen('slepian',lwin,TW,TC,nwin)

% Last modified by fjsimons-at-mit.edu, Dec 11th, 2000

switch funname
  case 'hermite'
    [lwin,R,nwin]=deal(varargin{:});
    [H,V,K,tax]=choice(R,lwin);
    if isempty(nwin); nwin=K; end
    stru={ 'Hermite functions'};
    [totabstf,fax,tax]=doit(H,V,nwin,tax);
    if ~nargout
      plotit(tax,H,nwin,V,fax,totabstf,R,stru)
    end

  case 'slepian'
    [lwin,TW,TC,nwin]=deal(varargin{:});
    [H,V,fw,fc]=slepian(lwin,TW,TC,nwin);    
    stru={ 'Slepian wavelets'};
    [totabstf,fax,tax]=doit(H,V,nwin,0:lwin-1);
    if ~nargout
      plotit(tax,H,nwin,V,fax,totabstf,[fw fc],stru)
    end
end  

% Ouptut if needed
posout={ 'totabstf','fax','tax'};
for index=1:nargout
  eval([ 'varargout{index}=',posout{index},';'])
end

% Supporting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [totabstf,fax,tax]=doit(H,V,nwin,tax)
for index=1:nwin
  [abstf(:,:,index),fax,tax]=wigner(H(:,index),indeks(diff(tax),1),tax);
  abstf(:,:,index)=abstf(:,:,index)*V(index);
end
totabstf=sum(abstf,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotit(tax,H,nwin,V,fax,totabstf,Rorfwfc,stru)
subplot(221)
plot(tax,H(:,1:nwin),'b-')
axis tight
hold on; plot(xlim,[0 0],'k'); hold off
title(stru,'FontS',15)

subplot(222)
plot(0:nwin-1,V(1:nwin),'-o','LineW',2)
title('Eigenvalues','FontS',15)
ylim([0 1.1])

subplot(212)
% Now it depends if it's Hermite or Slepian
if size(Rorfwfc)==1
  contourf(tax,2*pi*fax,totabstf);  
  axis image xy    
  ylim([0 (Rorfwfc+1)])
  title([ 'Concentration radius: ',sprintf('%i',Rorfwfc)],'FontS',15)
  ylabel('2\pi{f}')
else
  contourf(tax,fax,totabstf); hold on
  fw=Rorfwfc(1);  
  fc=Rorfwfc(2); 
  plot(xlim,[fc fc],'k--','LineW',2)
  plot([range(xlim)/2 range(xlim)/2],[fc+fw fc-fw],'k--','LineW',2)
  title([ 'Center frequency ',sprintf('%6.3f',fc),...
	'; Bandwidth ',sprintf('%6.3f',fw)],'FontS',15)
  ylabel('f')
  ylim([0 2*fc])
end

xlabel('Time')
colorbar('ver')





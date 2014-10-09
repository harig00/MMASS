function loading1
% LOADING1
%
% Plots standard admittance and coherence curves for various
% scenarios of surface and subsurface loading on two interfaces.
%
% Demo for program LOADING.
% 
% See also MCKENZIE

[xfaks,yfaks,fnx,fny,xsint,ysint]=...
    fftaxis([128 128],[128 128],[2000 2000]);
kx=2*pi*xfaks; ky=2*pi*yfaks; 

% Densities, depths and elastic thickness
defval('R',[0 2600 2900 3300])
defval('Te',35);
defval('T',[0 15 35]);

F2=[0 0.5 0.7 0.8 1];
F1=1-F2;
F3=zeros(size(F2));

for index=1:length(F2) 
  [jk1,jk2,Z(:,index),C2(:,index),k]=...
      loading(R,T,Te,[],kx,ky,...
      [F1(index) F2(index) F3(index)]);
end

xel=2*pi./[2001 100];
xtk={ '0.004' ' ' ' ' '0.007' ' ' ' ' '0.01' '0.02' '0.03' ' ' '0.05' ' '};
xtie=[0.004:0.001:0.01 0.02:0.01:0.06];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(1)=subplot(121);
sl(:,1)=semilogx(k,Z);
keyboard
set(gca,'Xdir','Norm','Ydir','Rev')
xlim(xel)
set(ah(1),'xtick',xtie,'XtickL',xtk)
ylim([-0.2 0.05])
xt(1)=xtraxis1d(ah(1));
xl(1)=xlabel('Bouguer Admittance');
axes(ah(1))
yl(1)=ylabel('Z_B(\lambda) (mgal/m)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(2)=subplot(122);
sl(:,2)=semilogx(k,C2);
set(gca,'Xdir','Norm')
xlim(xel)
set(ah(2),'xtick',xtie,'XtickL',xtk)
ylim([0 1.1])  
xt(2)=xtraxis1d(ah(2));
xl(2)=xlabel('Bouguer Coherence');
axes(ah(2))
yl(2)=ylabel('\gamma^2_B(\lambda)');

longticks(ah)
longticks(xt)

set(sl,'LineW',2)
set([xl yl],'FontS',15)

movev([ah xt],-.02)

shrink([ah xt],1,2)

axes(ah(1))
xlx(1)=xlabel('k=2\pi/\lambda (rad/km)','FontS',12);
lg(1)=legend(sprintf('T_e= %i km',Te),3);

axes(ah(2))
xlx(2)=xlabel('k=2\pi/\lambda (rad/km)','FontS',12);
lg(2)=legend(cellstr(eval(['[' sprintf('''F2= %2.1f'';',F2) ']'])),3);

fig2print(gcf,'landscape')
id

set(xt,'xgrid','on','ygrid','on')
axes(xt(1))
axes(xt(2))
axes(lg(1))
axes(lg(2))

figdisp

function varargout=loading2(Te,T)
% [Hs,Gb]=LOADING2(Te,T)
%
% Illustrates the use of the function LOADING
%
% Last modified by fjsimons-at-alum.mit.edu, May 13th 2003

defval('Te',25)
defval('T',[     0   15   35])
defval('R',[0 2600 2900 3300])
ddir= '/u/fjsimons/MyPapers/2003/JGR-2003/DATA/COHSYNTHETIC/';
load(fullfile(ddir,'Xmidcpy'))
bs=length(gravcpy);

% The larger this distance, the better it works, of course 
lenx=3000;
leny=3000;
x=linspace(0,lenx,size(gravcpy,2));
y=linspace(0,leny,size(gravcpy,1));
[xfaks,yfaks,fnx,fny,xsint,ysint]=fftaxis([bs bs],[bs bs],[leny lenx]);
% Definitely use this thing so loading comes up with REAL results
[K,kx,ky]=knum2(size(gravcpy),[leny lenx]);

% Create loading topography at three interfaces
AT(:,:,1)=scale(topocpy,[-0.5 0.5]);
% These are two uncorrelated processes, for sure;
% the data come from somewhere else
AT(:,:,2)=scale(gravcpy,[-3 3]);
%AT(:,:,2)=zeros(size(topocpy));
AT(:,:,3)=zeros(size(topocpy));

[Hs,Gb,Z,C2,k,ZK,C2K,f]=loading(R,T,Te,AT,kx,ky);

% Plot these things
clf
ah=krijetem(subnum(4,2));
fig2print(gcf,'landscape')

axes(ah(1))
imagesc(x,y,AT(:,:,1))
cb(1)=colorbar;
axis image
t(1)=title('Topography at z=0');

axes(ah(2))
imagesc(x,y,AT(:,:,2))
cb(2)=colorbar;
axis image
t(2)=title(sprintf('Topography at z= %2.0i km',T(2)));

axes(ah(3))
imagesc(x,y,Hs)
cb(3)=colorbar;
axis image
t(3)=title('Final Topography');

axes(ah(4))
imagesc(x,y,Gb)
cb(4)=colorbar;
axis image
t(4)=title('Bouguer Gravity Anomaly');

axes(ah(5))
imagesc(kx,ky,real(ZK))
cb(5)=colorbar;
axis image

axes(ah(6))
imagesc(kx,ky,C2K)
caxis([0 1])
cb(6)=colorbar;
axis image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xel=2*pi./[lenx 2*xsint];
resolutionlengthNW4=2*pi*4/lenx;
disp(sprintf('Resolution for NW=4 equals %8.3g in wavenumber',resolutionlengthNW4))
xtk={ '0.004' ' ' ' ' '0.007' ' ' ' ' '0.01' '0.02' '0.03' ' ' '0.05' ' '};
xtie=[0.004:0.001:0.01 0.02:0.01:0.06];

axes(ah(7))
sl(:,1)=semilogx(k,real(Z),'Color','k','LineW',2);
set(ah(7),'Xdir','Norm','Ydir','Rev')
ylim([-0.2 0.05])
xlim(xel)
set(ah(7),'xtick',xtie,'XtickL',xtk)
yl(1)=ylabel('Re Z_B(\lambda)');
xlx(1)=xlabel('k=2\pi/\lambda (rad/km)','FontS',12);
xt(1)=xtraxis1d(ah(7));
grid on
xl(1)=xlabel('Bouguer Admittance');
axes(ah(7))
lg(1)=legend(cellstr(eval(['[' sprintf('''f= %4.2f'';',f) ']'])),3);
axes(xt(1))
axes(lg(1))

axes(ah(8))
sl(:,2)=semilogx(k,C2,'Color','k','LineW',2);
set(gca,'Xdir','Norm')
ylim([0 1.1])  
xlim(xel)
set(ah(8),'xtick',xtie,'XtickL',xtk)
yl(2)=ylabel('\gamma^2_B(\lambda)');
xlx(2)=xlabel('k=2\pi/\lambda (rad/km)','FontS',12);
xt(2)=xtraxis1d(ah(8));
grid on
xl(2)=xlabel('Bouguer Coherence');
axes(ah(8))
lg(2)=legend(cellstr(eval(['[' sprintf('''Te= %3.1f'';',Te) ']'])),3);
axes(xt(2))
axes(lg(2))

set([t xl],'FontS',15)
set(sl,'LineW',2)
fig2print(gcf,'landscape')

nolabels(ah(1:4))
drawnow 

[FX,FY,SX,SY,SXY,coh20,E,W,coh2var,Z20]=mtm4(Hs,Gb,2);
[FX,FY,SX,SY,SXY,coh25,E,W,coh2var,Z25]=mtm4(Hs,Gb,2.5);
[FX,FY,SX,SY,SXY,coh30,E,W,coh2var,Z30]=mtm4(Hs,Gb,3);
[FX,FY,SX,SY,SXY,cohNW,E,W,coh2var,ZNW]=mtm4(Hs,Gb,4);

% Compare coherences
[knum,av20]=isav(coh20,kx,ky);
[knum,av25]=isav(coh25,kx,ky);
[knum,av30]=isav(coh30,kx,ky);
[knum,avNW]=isav(cohNW,kx,ky);

% Compare admittances
[knum,Zav20]=isav(Z20,kx,ky);
[knum,Zav25]=isav(Z25,kx,ky);
[knum,Zav30]=isav(Z30,kx,ky);
[knum,ZavNW]=isav(ZNW,kx,ky);

axes(ah(7))
hold on
pmZ(1)=semilogx(knum,scale(real(Zav20),minmax(real(Z))),'v-');
pmZ(2)=semilogx(knum,scale(real(Zav25),minmax(real(Z))),'^-');
pmZ(3)=semilogx(knum,scale(real(Zav30),minmax(real(Z))),'s-');
pmZ(4)=semilogx(knum,scale(real(ZavNW),minmax(real(Z))),'o-');
axes(xt(1))
axes(lg(1))

axes(ah(8))
hold on
pm(1)=semilogx(knum,av20,'v-');
pm(2)=semilogx(knum,av25,'^-');
pm(3)=semilogx(knum,av30,'s-');
pm(4)=semilogx(knum,avNW,'o-');

axes(xt(2))
axes(lg(2))

set(ah([5 6]),'Xaxisl','t')
moveh(cb,-.05)
movev([ah([5 6]) cb([5 6])],.02)

disp(sprintf('print -depsc %s/EPS/%s_%3.3i_%3.3i',...
	     getenv('GIF'),'loading2',Te,ceil(f*100)))

nargs={ 'Hs','Gb'};
for index=1:nargout
  varargout{index}=eval(nargs{index});
end



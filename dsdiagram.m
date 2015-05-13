function dsdiagram
% Geometry
%
% Simons & Dahlen (2007), Figure 1
%
% Last modified by fjsimons-at-alum.mit.edu, 03/28/2007

% Which vector to plot
ang=40;
% Down to this z level for the projection
lz=-0.2;
% Rotation of the geodesic
rotg=-30;
% Rotation of the X-axis
rots=10;

clf
[ah,ha]=krijetem(subnum(1,4));

% FIRST FIGURE PANEL ------------------------------
axes(ah(1))
[h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
eq(1)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
hold on
[h,cord]=circ(1,[0 3*pi/2]); delete(h)
eqd(2)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
[h,cord]=circ(1); delete(h)
ol(1)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));

% Plot the X-axis
[axt(1),axh(1)]=plotx;

% Plot the Y-axis
[axt(2),axh(2)]=ploty;

% Plot the Z-axis
[axt(3),axh(3)]=plotz;

% Plot the random vector
ft=[1 1]; 
vax=arrow(0,0,0,0.85,'v',2,ang,ft);
xdv=get(vax,'Xdata'); ydv=get(vax,'Ydata'); delete(vax)
axt(4)=plot3(zeros(size(ydv{1})),xdv{1}-ft(1),ydv{1}-ft(2));
%axh(4)=plot3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2));
ydv{2}(4)=ydv{2}(1);
xdv{2}(4)=xdv{2}(1);
axh(4)=fill3(zeros(size(ydv{2})),xdv{2}-ft(1),ydv{2}-ft(2),'k');
% Plot first arclength
[h,cord]=circ(0.5,[pi/2-ang/180*pi pi/2]); delete(h)
arcl(1)=plot3(zeros(size(cord(:,1))),cord(:,1),cord(:,2));
% Plot the arrow head
fudg=10; % Used to be 4
parh=plot3(0,cord(1+fudg,1),cord(1+fudg,2),'k^');

% Plot projections
proj(1)=plot3([0 0],repmat(max(xdv{1})-ft(1),1,2),...
	      [max(ydv{1})-ft(1) lz]);
proj(2)=plot3([0 0],[0 max(xdv{1})-ft(1)],[0 lz]);

% Plot second arclength
[h,cord]=circ(0.4,[-rots*pi/180 ang/180*pi]); delete(h)
arcl(2)=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
% Plot the second arrow head
parh(2)=plot3(cord(end,1)+0.025,cord(end,2)-0.05,0,'kv');
set(parh,'markerf','k','markere','k','markers',3)

onax=axis; onax=onax.*1.1;
viewpars(onax)

% SECOND FIGURE PANEL
axes(ah(2))
% Plot Earth outline
[h,cord]=circ(1); delete(h)
ol(3)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
hold on
% First blob
[lo,la]=blob; 
lo=scale(lo,(-60+[0 50])*pi/180);
la=scale(la,[40 70]*pi/180);
[X,Y,Z]=sph2cart(lo,la,1);
arcl(3)=fill3(X,Y,Z,grey);
% Second blob
[lo,la]=blob; 
lo=scale(lo,(30+[-40 0])*pi/180);
la=scale(la,[20 -20]*pi/180);
[X,Y,Z]=sph2cart(lo,la,1);
arcl(6)=fill3(X,Y,Z,grey);
curs=0;
if curs==1
  [phint,thp,php]=phicurve([pi/2-la(:) lo(:)],...
			   linspace(pi/2-max(la),pi/2-min(la),10));
  % Now get the great circle coordinates
  thp=thp'; php=php';
  for index=1:size(php,1)
    lola=grcircle(...
	[php(index,1) pi/2-thp(index,1)],...
	[php(index,2) pi/2-thp(index,2)],100);
    [X,Y,Z]=sph2cart(lola(:,1),lola(:,2),1);
    hasj(index)=plot3(X,Y,Z,'k');
  end
end
viewpars(onax)

% THIRD FIGURE PANEL
axes(ah(3))
% Plot Earth outline
[h,cord]=circ(1); delete(h)
ol(4)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
hold on
% Plot polar cap
ang1=30;
[h,cord]=circ(cos(ang1*pi/180),[-pi/2 pi/2]); delete(h)
eq(4)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang1*pi/180));
[h,cord2]=circ(1,[ang1 180-ang1]*pi/180); delete(h)
eqx=fill3([cord(:,1) ;  ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ;  ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang1*pi/180) ; ... 
	     ; cord2(:,2)]',...
	    grey);
viewpars(onax)

% FOURTH FIGURE PANEL --------------------------------------------
axes(ah(4))
[h,cord]=circ(1); delete(h)
ol(4)=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
hold on
viewpars(onax)
% Plot polar cap on TOP
ang2=20;
[h,cord]=circ(cos(ang2*pi/180),[-pi/2 pi/2]); delete(h)
eq(6)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang2*pi/180));
[h,cord2]=circ(1,[ang2 180-ang2]*pi/180); delete(h)
eqx=fill3([cord(:,1) ;  ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ;  ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang2*pi/180) ; ... 
	     ; cord2(:,2)]',...
	    grey);
% Plot bottom polar cap on TOP
ang2=-ang2;
% Make an equatorial half circle
[h,cord]=circ(cos(ang2*pi/180),[-pi/2 pi/2]); delete(h)
% But plot it at the right height
eq(3)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang2*pi/180));
% Define another boundary
[h,cord2]=circ(1,[-ang2 180+ang2]*pi/180); delete(h)
% Make this grey
eqx=fill3([cord(:,1) ; zeros(size(cord2(:,1)))]',...
	  [cord(:,2) ; cord2(:,1)]',...
	  [ones(size(cord(:,1)))*sin(ang2*pi/180) ; -cord2(:,2)]',...
	  grey);

viewpars(onax)

% Cosmetics & Annotations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
serre(ah(1:4),-1/3,'across')
%Next couple of lines for krijetem(subnum(2,2))
%serre(ah(1:2),1,'across')
%serre(ah(3:4),1,'across')
%serre(ha(1:2),0.6,'down')
%serre(ha(3:4),0.6,'down')
%set(ah,'camerav',5.5) 
set(ah,'camerav',4.75) 
% Next line not for 2,2
set(ah,'xlim',[-1.25 1.25],'ylim',[-1.25 1.25],'zlim',[-1.25 1.25])
set(findobj('Color','b'),'LineW',1,'Color','k')
set([eq(~~eq) proj(~~proj)],'LineS',':')

% Plot labels
axes(ah(1))
lpost=[1.9  -0.4    0  ;
      0.5    1.2    0  ;
      0      0.125  1.2;
      0      0.5    0.7];
axpo=[1 1 1 1];
ltxt={'x','y','z','\bfr'};
for i=1:length(axpo)
  axes(ah(axpo(i)))
  l(i)=text(lpost(i,1),lpost(i,2),lpost(i,3),ltxt{i});
end

axes(ah(1))
l(8)=text(0,0.2,0.6,'\theta');
l(9)=text(0.65,0.2,0.02,'\phi');
l(7)=text(0,0.4,-0.7,'\Omega');
axes(ah(2))
l(11)=text(0,0.53,-0.45,'R_2');
l(13)=text(0,-0.05,0.6,'R_1');
axes(ah(3))
l(14)=text(0,0.93,0.55,sprintf('%s','\Theta'));
axes(ah(4))
l(15)=text(0,1.0,0.4,sprintf('%s','\Theta'));
l(16)=text(0,1.0,-0.4,sprintf('%s-%s','\pi','\Theta'));

fig2print(gcf,'portrait')
figdisp([],[],[],0)

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewpars(onix)
defval('onix',[-2 2 -2 2])
% Set viewing parameters
% Really would need another rotation around x
view(90,17.5); axis equal
axis(onix); axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h]=plotx
% Plots x-axis for the unit sphere
ft=[1 1]; % To get it on the right plane
% Rotation of the X-axis
rots=10;
xax=arrow(0,0,1.75,0,'h',2,rots,ft);
xdx=get(xax,'Xdata'); ydx=get(xax,'Ydata'); delete(xax)
t(1)=plot3(xdx{1}-ft(1),ydx{1}-ft(2),zeros(size(ydx{1})));
% Regular arrow head
% h(1)=plot3(xdx{2}-ft(1),ydx{2}-ft(2),zeros(size(ydx{2})));
% Filled arrow head, but unrotated
ydx{2}(4)=ydx{2}(1);
xdx{2}(4)=xdx{2}(1);
h(1)=fill3(xdx{2}-ft(1),ydx{2}-ft(2),zeros(size(ydx{2})),'k');
% Filled and properly rotated arrow head
% rots1=rotz(-rots*pi/180);
% rots2=rotx(pi/2);
% rots3=rotz(rots*pi/180);
% done=(rots3*rots2*rots1*[xdx{2}-ft(1),ydx{2}-ft(2),zeros(size(ydx{2}))]')';
% h(1)=fill3(done(:,1),done(:,2),done(:,3),'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h]=ploty
% Plots y-axis for the unit sphere
yax=arrow(0,0,0,1.4,'v',3);
xdy=get(yax,'Xdata'); ydy=get(yax,'Ydata'); delete(yax)
t(1)=plot3(xdy{1},ydy{1},zeros(size(ydy{1})));
% h(1)=plot3(zeros(size(ydy{2})),ydy{2},xdy{2});
ydy{2}(4)=ydy{2}(1);
xdy{2}(4)=xdy{2}(1);
h(1)=fill3(zeros(size(ydy{2})),ydy{2},xdy{2},'k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h]=plotz
% Plots z-axis for the unit sphere
zax=arrow(0,0,0,1.4,'v',3);
xdz=get(zax,'Xdata'); ydz=get(zax,'Ydata'); delete(zax)
t(1)=plot3(zeros(size(ydz{1})),xdz{1},ydz{1});
% h(1)=plot3(zeros(size(ydz{2})),xdz{2},ydz{2});
ydz{2}(4)=ydz{2}(1);
xdz{2}(4)=xdz{2}(1);
h(1)=fill3(zeros(size(ydz{2})),xdz{2},ydz{2},'k');




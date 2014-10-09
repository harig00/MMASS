function loris1(N,eo,lon,lat,sc,fax)
% LORIS1(N,eo,lon,lat,sc,fax)
%
% INPUT:
%
% N        The power of the dyadic subdivision [defaulted]
% eo       0 even number of points [default]
%          1 odd number of points
% lon,lat  The location of the viewing platform [defaulted]
% sc       0 regular cubed sphere [default]
%          1 superchunk cubed sphere
% fax      Axis scaling [defaulted]
%
% Makes a cute plot of the cubed sphere in three-dimensional space.
%
% Last modified by fjsimons-at-alum.mit.edu, 08/30/2011

% Set the defaults
defval('N',4)
% For GJI this was lon=-40 lat=65
%defval('lon',-40)
%defval('lat',65)
% For SPIE this was lon=30 lat=-5
defval('lon',30)
defval('lat',-5)
defval('actprint',0)
defval('eo',0)
defval('sc',0)
% For SPIE this was 1
fax=1;
% For GJI this was 0.8
defval('fax',0.8);

% Create the cubed sphere
Nor=N;
[x,y,z,J,N]=cube2sphere(N,[],[],[],eo,0);

if sc==1
  % Create the superchunk sphere and reassign N
  [x2,y2,z2,J,N2]=cube2sphere(Nor*2,[],[],[],eo,1);
end

clf
% Set view angles ahead of time as an explicit longitude and latitude
[xv,yv,zv]=sph2cart(lon*pi/180,lat*pi/180,1);

% Plot the visible continents, which depends on the view angle
[a,h,XYZ]=plotcont([0 90],[360 -90],3); delete(h); 
% Inner product selectivity
yes=[xv yv zv]*XYZ'>0; XYZ=XYZ(yes,1:3);
% This protection from jumps is straight from PLOTCONT
xx=XYZ(:,1); yy=XYZ(:,2); zz=XYZ(:,3);
d=sqrt((xx(2:end)-xx(1:end-1)).^2+(yy(2:end)-yy(1:end-1)).^2);
dlev=3; pp=find(d>dlev*nanmean(d));
nx=insert(xx,NaN,pp+1); ny=insert(yy,NaN,pp+1); nz=insert(zz,NaN,pp+1); 
% And then finally do it
skl=0.99;
pc=plot3(nx(:)*skl,ny(:)*skl,nz(:)*skl,'k-');
hold on

% Plot a single "main" face, all visible
% Don't do mesh as this is not see-through
% pm=mesh(x(:,:,1),y(:,:,1),z(:,:,1),ones(N,N),'edgecolor','k');
if sc==0
  % Figure out in what panel lies the center view point
  [~,~,mf]=sphere2cube(lon,lat);
  horz=[x(:,:,mf)  y(:,:,mf)  z(:,:,mf) ];
  vert=[x(:,:,mf)' y(:,:,mf)' z(:,:,mf)'];
  pm=plot3(horz(:,1:N),horz(:,N+1:2*N),horz(:,2*N+1:3*N),'k');
  um=plot3(vert(:,1:N),vert(:,N+1:2*N),vert(:,2*N+1:3*N),'k');
end

% Then plot the ribs only
for ind=1:6
  % These will be the lines plotted, with some redundancy
  horz=[x(:,[1 N],ind)  y(:,[1 N],ind)  z(:,[1 N],ind) ];
  vert=[x([1 N],:,ind)' y([1 N],:,ind)' z([1 N],:,ind)'];
  % But restrict them to the viewable area
  horz(repmat(reshape([xv yv zv]*[reshape(horz,2*N,[])]'<0,N,2),1,3))=NaN;
  vert(repmat(reshape([xv yv zv]*[reshape(vert,2*N,[])]'<0,N,2),1,3))=NaN; 
  % Then plot them for real
  p{ind}=plot3(horz(:,1:2),horz(:,3:4),horz(:,5:6),'k');
  u{ind}=plot3(vert(:,1:2),vert(:,3:4),vert(:,5:6),'k');
end

if sc==1
  % Now hold on and plot the superchunk ribs also
  x=x2; y=y2; z=z2; N=N2;
  for ind=1:6
    % These will be the lines plotted, with some redundancy
    horz=[x(:,[1 N],ind)  y(:,[1 N],ind)  z(:,[1 N],ind) ];
    vert=[x([1 N],:,ind)' y([1 N],:,ind)' z([1 N],:,ind)'];
    % But restrict them to the viewable area
    horz(repmat(reshape([xv yv zv]*[reshape(horz,2*N,[])]'<0,N,2),1,3))=NaN;
    vert(repmat(reshape([xv yv zv]*[reshape(vert,2*N,[])]'<0,N,2),1,3))=NaN; 
    % Then plot them for real
    p2{ind}=plot3(horz(:,1:2),horz(:,3:4),horz(:,5:6),'k');
    u2{ind}=plot3(vert(:,1:2),vert(:,3:4),vert(:,5:6),'k');
  end
end

% Now set (and verify the syntax) of the view 
view([xv,yv,zv]); [AZ,EL]=view;
disp(sprintf('Azimuth: %i ; Elevation: %i',round(AZ),round(EL)))

% And now plot an entire equatorial circle also
[xe,ye,ze]=sph2cart(linspace(0,2*pi,100),0,1);
xyze=[rotz(-lon*pi/180)*roty(-[90-lat]*pi/180)*[xe ; ye ; repmat(ze,1,length(ye))]]';
peq=plot3(xyze(:,1),xyze(:,2),xyze(:,3),'k');

% Cosmetics
% Change all of their colors
set(pc,'Color',grey,'LineW',1)
set(cat(1,p{:}),'Color','k','LineW',1)
set(cat(1,u{:}),'Color','k','LineW',1)
if sc==1
  set(cat(1,p2{:}),'Color','b','LineW',0.5)
  set(cat(1,u2{:}),'Color','b','LineW',0.5)
  % Save money for SPIE
  set(cat(1,p2{:}),'Color',grey,'LineW',0.5,'LineS','-')
  set(cat(1,u2{:}),'Color',grey,'LineW',0.5,'LineS','-')
end
% set(pm,'EdgeColor','k','LineW',0.5)
if sc==0
  set(pm,'Color','k','LineW',0.5)
  set(um,'Color','k','LineW',0.5)
end
set(peq,'Color','k','LineW',1)

fig2print(gcf,'portrait')

axis equal; axis([-1 1 -1 1 -1 1]*fax);
xl=xlabel('x'); yl=ylabel('y'); zl=zlabel('z');
set(gca,'xtick',[-fax 0 fax],'ytick',[-fax 0 fax],'ztick',[-fax 0 fax])
moveh(xl,-0.55); movev(xl,0.15)
moveh(yl,-0.2); movev(yl,0.2)
moveh(zl,-0.1); movev(zl,0.1)
axis off
% As if the next line did anything - it doesn't
bottom(pc,gca)

% Where is the North Pole?
hold on
pnp=plot3(0,0,1,'MarkerF','k','MarkerE','k','Marker','o');
delete(pnp)
hold off
% Where is the Viewing Axis?
hold on
pnp=plot3(xv,yv,zv,'MarkerF','k','MarkerE','k','Marker','o');
delete(pnp)
hold off

% Print it out
figna=figdisp([],sc,[],actprint);
if actprint==1
  system(sprintf('epstopdf %s.eps',figna));
end


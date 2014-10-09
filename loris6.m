function loris6(wav)
% LORIS6(wav)
%
% This is supposed to plot some wavelets on the cubed sphere
%
% INPUT:
%
% wav    The chosen construction, 'D2', 'D4' or 'D6', or 'CDF4.2'
%
% Last modified by fjsimons-at-alum.mit.edu, 01/30/2011

% Here you want to use the same values as in LORIS3
defval('N',7)
defval('wav','D4')
defval('J',N-3-strcmp(wav,'D6'))

% Plot
clf
[ah,ha,H]=krijetem(subnum(2,5));

% First the simple grid and the centers of squares
axes(ah(1))
[ph,W]=fridplotw(N,J); hold on
axis image off

% Sort somehow
[~,j]=sort(sqrt(sum(W.^2,1)));
W=W(:,j);

% Make sure we know what we'll get
vwlev1=cube2scale(N,[J J],1,1);

% Get the continents, short from PLOTCONT
[xic,etac]=deal([0,0]);%getcont(3,N);

% Now plot the wavelets
for index=2:length(ah)
  zskel=vwlev1(W(1,index-1),W(2,index-1));
  if strcmp(wav,'CDF4.2')
    delso=20/2^(J-zskel);
  elseif strcmp(wav,'D4')
    delso=15/2^(J-zskel);
  end
  
  ttt(index-1)=plotit(ah(index),W(1,index-1)+3,W(2,index-1)+3,N,J,wav,delso);
  movev(ttt(index-1),-10)
  hold on
  p(index-1)=plot(etac,xic,'k');
  hold off
  t(index-1)=text(77,10,sprintf('scale %i',zskel));
  tt(index-1)=text(5,9,sprintf('N = %i',N));
  axis image 
  axes(ah(1))
  % Remember which is first here
  % slight upping
  if index<=5
    extr=0.75;
  else
    extr=0;
  end
  l(index-1)=text(W(2,index-1)+3,W(1,index-1)+3+extr,letter(index-1));
end
hold off

% Cosmetics
set(l,'FontS',4,'HorizontalA','c','VerticalA','m')
set([t ttt],'FontS',8)
set([tt],'FontS',6)
[bh,th]=label(ah(2:end),'ur',8);
noticks(ah(2:end))
nolabels(ah(2:end))
movev(H(1,:),-.275)
fig2print(gcf,'portrait')
figdisp([],wav,[],0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ttit=plotit(aha,w1,w2,N,J,wav,delso)

defval('wav','D4')

% First the zeros
x=zeros(2^N,2^N); 

% Then the one one
x(w1,w2)=1;

defval('xver',0)

switch wav
  case 'D2'
   % Calculate this thing
   f=angularD2WT(x,[J J],'inverse',1);
   
   if xver==1
     % Do I get the same thing back?
     xx=angularD2WT(f,[J J],'forward',1);
   end
  case 'D4'
   % Calculate this thing
   f=angularD4WT(x,[J J],[0 0],'inverse',1);
   
   if xver==1
     % Do I get the same thing back?
     xx=angularD4WT(f,[J J],[0 0],'forward',1);
   end
  case 'D6'
   % Calculate this thing
   f=angularD6WT(x,[J J],[0 0],'inverse',1);
   
   if xver==1
     % Do I get the same thing back?
     xx=angularD6WT(f,[J J],[0 0],'forward',1);
   end
 case 'CDF4.2'
  % I temporarily rely on Ignace's calculation here
  load /u/fjsimons/MyPapers/2011/GJI-2011/cdf42examplesafrica.mat
  indices=[96 48 24 12  4 32 16  8  4 96 48 24 12 ;
           96 48 24 12  4 96 48 24 12 32 16  8  4];
  itsit=find(indices(1,:)==w1 & indices(2,:)==w2);
  f=africanwavelets(:,:,itsit);
end

try
  % Check the norm of the reconstruction
  difer(x-xx,[],[],NaN)
end

defval('xver',0)
if xver==1
  % Check orthgonality by calculating any other one and computing the inner
  % product. Here check the orthonormality and one vanishing moment; they
  % should have two each.
  disp(sprintf('Orthohonormality %g',sum(sum(f.^2))))
  disp(sprintf('Vanishing moment 0 %g',sum(sum(f))))
  xi=linspace(-pi/4,pi/4,size(f,2));
  eta=linspace(-pi/4,pi/4,size(f,1));
  [XI,ETA]=meshgrid(xi,eta);
  disp(sprintf('Vanishing moment 1 %g',sum(sum(XI.*f))))
  disp(sprintf('Vanishing moment 1 %g',sum(sum(ETA.*f))))
  disp(sprintf('Vanishing moment 1 %g',sum(sum(XI.*ETA.*f))))
  disp(sprintf('Vanishing moment 2 %g',sum(sum(XI.*XI.*f))))
  disp(sprintf('Vanishing moment 2 %g',sum(sum(ETA.*ETA.*f))))
end

% Plot it
axes(aha)
imagefnan([1 1],[2^N 2^N],f,[],halverange(f,75))

if xver==1
  clf
  ff=zeros([size(f) 6]);
  plotoncube(ff,'2D')
  ff(:,:,3)=f;
  plotoncube(ff,'2D')
  hold on
  plotcont([],[],9)  
end

% Probably time to calculate the dominant degree in the spherical
% harmonic expansion also?
hold off
witsj=3;
% This is by interpolation to a regular grid; you might want to avoid
% this and go straight to the spherical harmonics
[lons,lats,vd,lon,lat]=cube2lola(f,witsj); vd(isnan(vd))=0;
% Figure out the dominant wavelength
%sdl,l]=plm2spec(xyz2plm(vd));
%semilogy(l,sdl)
%a,b]=max(sdl);
% This doesn't turn out to be too meaningful
%tit=title(sprintf('max l = %i',l(b),);
% And what is the epicentral distance of the patch, instead
% Remember that elsewhere everything is flipped!
[i,j]=find(~~flipud(f));
% Find the x,y coordinates of corners and the center
cone=[lon(min(i),min(j)) lat(min(i),min(j))]*180/pi;
ctwo=[lon(min(i),max(j)) lat(min(i),max(j))]*180/pi;
cthr=[lon(max(i),min(j)) lat(max(i),min(j))]*180/pi;
cfor=[lon(max(i),max(j)) lat(max(i),max(j))]*180/pi;
cen=[lon(min(i)+round(diff(minmax(i))/2),...
	 min(j)+round(diff(minmax(j))/2)) ...
     lat(min(i)+round(diff(minmax(i))/2),...
	 min(j)+round(diff(minmax(j))/2))]*180/pi;
% Epicentral distances of the corners to the center
[gcdkm,dlto]=grcdist(cen,[cone ; ctwo ; cthr ; cfor]);
defval('delso',mean(dlto)/2);
% Plot a first circle
[lon2,lat2]=caploc(cen,delso);
[xic,etc]=lola2face(lon2,lat2,N,witsj);
hold on
p(1)=plot(etc,xic,'b-');
% Plot a second circle
[lon2,lat2]=caploc(cen,delso/2);
[xic,etc]=lola2face(lon2,lat2,N,witsj);
p(2)=plot(etc,xic,'r-');

set(p(1),'linew',0.25,'Color','k','LineS','-')
set(p(2),'linew',0.25,'Color','k','LineS','-')

if xver==1
  % Plot the center also - note that we have a xi-eta switch
  [xic,etc]=lola2face(cen(1),cen(2),N,witsj);
  p(3)=plot(etc,xic,'ro');
  % Plot the four corners also
  [xic,etc]=lola2face(cone(1),cone(2),N,witsj);
  p(4)=plot(etc,xic,'go');
  [xic,etc]=lola2face(ctwo(1),ctwo(2),N,witsj);
  p(5)=plot(etc,xic,'go');
  [xic,etc]=lola2face(cthr(1),cthr(2),N,witsj);
  p(6)=plot(etc,xic,'go');
  [xc,etc]=lola2face(cfor(1),cfor(2),N,witsj);
  p(7)=plot(etc,xic,'go');
end

ttit=title(sprintf('%s = %g%s, %g%s','\Delta',...
		   round(delso/2*10)/10,str2mat(176),...
		   round(delso*10)/10,str2mat(176)));

if xver==1
  % Check this out on the map!
  clf
  plotplm(vd,[],[],4)
  hold on
  % plot(lon*180/pi,lat*180/pi,'k.'); 
  [~,pc]=plotcont; hold on; 
  pu=plot(lon(i,j)*180/pi,lat(i,j)*180/pi,'w.'); axis([0 80 -60 60])
  pp(1)=plot(cone(1),cone(2),'bo');
  pp(2)=plot(ctwo(1),ctwo(2),'bo');
  pp(3)=plot(cthr(1),cthr(2),'bo');
  pp(4)=plot(cfor(1),cfor(2),'bo');
  pp(5)=plot(cen(1),cen(2),'ro');
  set(pp,'MarkerF','y')
  set(pc,'LineW',2)
  % Try this now also, NOPE, this is too heavy
  % lon=lon(:,:,3)*180/pi;
  % lat=lat(:,:,3)*180/pi;
  % [sdl2,l2]=plm2spec(xyz2plm(f(:),128,[],lat(:),lon(:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xic,etc]=getcont(in,N)
% Should try this in LORIS3 also
[~,h,XYZ]=plotcont([],[],9); delete(h)
% in is the face index
xic =[XYZ(1,in)]; etc=[XYZ(2,in)];
% This isn't precise yet, should try an overlay later
xic=scale([xic{:} -pi/4 pi/4],[1 2^N]); xic=xic(1:end-2);
etc=scale([etc{:} -pi/4 pi/4],[1 2^N]); etc=etc(1:end-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xic,etc]=lola2face(lon,lat,N,witsj)
[xi,eta,fid]=sphere2cube(lon,lat);
difer(fid-witsj,[],[],NaN)
xic=scale([xi{witsj} ; -pi/4 ; pi/4],[1 2^N]); xic=xic(1:end-2);
etc=scale([eta{witsj}; -pi/4 ; pi/4],[1 2^N]); etc=etc(1:end-2);

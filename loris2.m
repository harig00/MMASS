function loris2(N,L,colmap,meth,fs)
% LORIS2(N,L,colmap,meth,fs)
%
% Makes a cute plot of the cubed sphere in the plane.
%
% INPUT:
%
% N        The power of the dyadic subdivision [defaulted]
% L        Spherical harmonic degree limit [defaulted]
% colmap   Colormap string {'demmap' | 'kelicol' | 'sergeicol'}
% meth     The option called in PLOTONCUBE [defaulted]
% fs       Fonts size of the color bar in pixels [defaulted]
%
% Last modified by fjsimons-at-alum.mit.edu, 08/29/2010

% Set the defaults
defval('N',7) % 8 is now possible too
defval('L',ceil(2^(N+1)))
defval('colmap','sergeicol');
defval('meth',1)
defval('fs',8)

clf

% Load or make the data
fname=fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS',...
	       sprintf('loris2_%i_%i.mat',N,L));

if exist(fname,'file')==2
  load(fname)
else
  % Load Earth's topography
  lmcosi=rindeks(fralmanac('GTM3AR','SHM'),1:addmup(L));
  
  % Perform the transform with standard inputs
  v=plm2cube(lmcosi,N);

  % Save
  save(fname,'v','N','L')
end

clf
% Plot in the two-dimensional plane - no reordering needed
[dem,dax,ziro]=eval(colmap);

% Hardwire so that the zero crossing is visually pleasing
dax=[-7473 5731].*[1.025 0.975];

% Make the plot
ah=gca;
[p,pg,txtyw]=plotoncube(v,'2D',meth,[],[],[],dax,colmap,0);

if meth~=3
  % Plot the continents on top
  [a,pc]=plotcont([],[],9);
  set(pc,'lines','-','marker','none')
  % Put some diagnostics where they belong
  xloc=-0.75;
  t(1)=text(xloc,1.25,sprintf('N = %i',N));
  t(2)=text(xloc,1.00,sprintf('L = %i',L));
  xloc=4.05;
  t(3)=text(xloc,1.50,sprintf('max = %5i m',round(max(v(:)))));
  t(4)=text(xloc,1.25,sprintf('med = %5i m',round(median(v(:)))));
  t(5)=text(xloc,1.00,sprintf('min = %5i m',round(min(v(:)))));
end

% Labels and legends
warning off MATLAB:quiver:DeprecatedV6Argument
xnum=12.5;
ynum=10.5;
lnum=10;
pv=[5 4 6 2 3 1];
ct={'x+','z-','y+','x-','z+','y-'};
for in=1:6
  xwid=txtyw(in,3);
  ywid=txtyw(in,4);
  [bh,th]=boxtex(txtyw(in,1:2)+[xwid/xnum ywid/-ynum],...
		 ah,pv(in),fs+1,2);
  [bhh,thh]=boxtex(txtyw(in,1:2)+[xwid-xwid/xnum ywid/-ynum],...
		 ah,'z+',fs+1,2,[],1.2);
  % If not this way the boxes are unequal in size
  set(thh,'string',ct{in})
  arx=txtyw(in,1)+xwid/2/xnum;
  ary=txtyw(in,2)-ywid+ywid/2/ynum;
  tx{in}=arrow(arx,ary,0,ywid/lnum,'v');
  ty{in}=arrow(arx,ary,xwid/lnum,0,'h');
  ttx(in)=text(arx+xwid/lnum+xwid/lnum/10,ary+xwid/lnum/10,'\xi')';
  tty(in)=text(arx,ary+ywid/lnum-ywid/lnum/10,'\eta');
end
% Red will stand out
set(cat(1,tx{:}),'color','r')
set(cat(1,ty{:}),'color','r')
set([ttx tty],'color','r')
set(tty,'VerticalA','bottom')
set(ttx,'HorizontalA','left')
set(ttx,'VerticalA','mid')
set([ttx tty],'FontS',fs)
warning on MATLAB:quiver:DeprecatedV6Argument

% Add a nice colorbar
colpos=[0.5616    0.1714+0.025    0.3143    0.0298];
[cb,xcb]=addcb(colpos,dax,dax,colmap,range(dax)/5,0);
set(cb,'fonts',fs)
set(xcb,'string','topography (m)','fonts',fs)
% Be very explicit
set(cb,'xtick',[-7500:3750:3750])
set(cb,'xtickl',round(get(cb,'xtick')))

% Leave with the old axis active
axes(ah)

% Pure cosmetics
fig2print(gcf,'portrait')

% Print it out
figdisp([],[],[],0)

% Special for SPIE2011
delete(t)
delete(cb)
delete(cat(1,p{:}))
set(cat(1,tx{:}),'color','k')
set(cat(1,ty{:}),'color','k')
set([ttx tty],'color','k')

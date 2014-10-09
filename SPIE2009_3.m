function SPIE2009_3(agu,L,J)
% SPIE2009, Figure 1
%
% A plot of the L=72 circular Slepian functions concentrated on the
% Bangui anomaly or the Sumatran earthquake.
% 
% INPUT
%
% agu    1 Adaptation for AGU (default: 0)
%        2 Adaptation for ALUMNI
% L      Bandwidth (default: 72)
% J      Numbers shown (default: the Shannon number, or minimum 30)
%
% SEE ALSO:
%
% SPIE2009_8
%
% Last modified by fjsimons-at-alum.mit.edu, 05/26/2011

defval('agu',0)
defval('lox','ll')
defval('fozo',6)

% Geometry of the Slepian window
defval('L',72);
TH=18;
phi0=18;
theta0=85;
omega0=0;
defval('J',max(30,round(spharea(TH,1)*(L+1)^2)));
degres=1/4;
ofs=10;

if agu==2
  TH=10;
  theta0=85;
  phi0=95;
  omega0=0;
  % Calculate the circle with ofs degrees offset
  ofs=0;
  degres=1/4;
end

% Calculate the circle with 10 degrees offset to hit Africa right
[lon2,lat2]=caploc([phi0+ofs 90-theta0],TH);

% Watch the halving of the L which is a bit of an artificial fix
% TO WHAT? I CAN NO LONGER REMEMBER
fname=fullfile(getenv('IFILES'),'GLMALPHAPTO',...
	       sprintf('SPIE-%i-%i-%i-%i-%i.mat',...
		       TH,L/2,phi0,theta0,omega0));

% Define the grid
theta=[theta0-2*TH:degres:theta0+2*TH];
phi=[phi0-4*TH:degres:phi0+4*TH];
c11cmn=[min(phi) 90-(theta0-2*TH) max(phi) 90-(theta0+2*TH)];

if exist(fname)==2
  load(fname)
else
  [Gar,V]=galphapto(TH,L,phi0,theta0,omega0,theta/180*pi,phi/180*pi,J);
  save(fname,'Gar','V')
end

% Seems like by now it should have done this automatically, must be a fix
V=sort(V,'descend');

% Trick for the continents to show up nice
[dr,lola]=maprotate([],[phi0-2*TH 90 360+(phi0-2*TH) -90]);

% Make the plot
clf
[ah,ha,H]=krijetem(subnum(6,5));
for index=1:length(ah)-[agu~=0]*15
  axes(ah(index))
  drexp=reshape(Gar(index,:),length(theta),length(phi));
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(drexp),...
	    kelicol,halverange(drexp));
  hold on; plot(lola(:,1),lola(:,2),'k'); 
  plot(lon2-ofs,lat2,'k')
  hold off
  axis([phi0-2*TH phi0+2*TH 90-(theta0+2*TH) 90-(theta0-2*TH)])

  if agu~=0
    t{index}=sprintf('%s =%7.3f','\lambda',V(index));
    [bh(index),th(index)]=boxtex(lox,ah(index),t{index},fozo,[],[],1.1);
  else
    tl(index)=title(sprintf('%s_{%i} = %9.6f','\lambda',index,V(index)));
  end
  drawnow
end

% Cosmetics (best after KEYBOARD or running the whole thing once)
set(ah,'FontS',8)
longticks(ah,1/2)
fig2print(gcf,'tall')

keyboard

switch agu
 case 1
  set(ah,'xtick',[0:20:40],'ytick',[-20:20:40])
  deggies(ah)
  serre(H',1/2,'down')
  serre(H,1/2,'across')
  nolabels(ha(7:end),2)
  nolabels(ah(1:10),1)
  delete(ah(16:end))
  figdisp('agu2009_1',[],[],0)
 case 2
  set(ah,'xtick',[80:15:110],'ytick',[-10:10:20])
  deggies(ah)
  serre(H',1/2,'down')
  serre(H,1/2,'across')
  nolabels(ha(7:end),2)
  nolabels(ah(1:10),1)
  delete(ah(16:end))
  figdisp('alumni2011_1',[],[],0)
 case 0
  movev(tl,-3)
  set(tl,'FontS',8)
  set(ah,'xtick',[0:20:40],'ytick',[-20:20:40])
  deggies(ah)
  serre(H',1/3,'down')
  nolabels(ha(7:end),2)
  nolabels(ah(1:end-5),1)
  figdisp([],[],[],0)
end

function SPIE2009_4
% SPIE2009_4
%
% A plot of the L=60 expanded gravitational potential perturbations due
% to some large end-member yearthquakes.
%
% Last modified by fjsimons-at-alum.mit.edu, 06/03/2010

disp('Until normalization is fixed we will not have a good plot')

% Produce fictitious earthquake at the equator
defval('deplonlat',[31 95 5])
defval('degres',0.5)
defval('Lup',60)

% Make the moment tensors, see Dahlen & Tromp Table 5.1
M{1}=[ 1  0 -1  0 0 0 ]/sqrt(2); 
M{2}=[ 0  0  0  0 1 0 ]/sqrt(2); 
M{3}=[ 0  1 -1  0 0 0 ]/sqrt(2); 
M{4}=[ 1 -1  0  0 0 0 ]/sqrt(2); 
M{5}=[ 0  0  0  1 0 0 ]/sqrt(2); 
M{6}=[ 0  0  0  0 0 1 ]/-sqrt(2); 

% Produce the first of the potential perturbations, collect eigenfunctions
[p{1},pr{1},Sstrain]=eqpotential(M{1},deplonlat,Lup);

% Produce the rest, reuse eigenfunctions
for ind=2:length(M)
  [p{ind},pr{ind},Sstrain]=eqpotential(M{ind},deplonlat,Lup,Sstrain);
end

% Then do the plotting, for speed only expand locally
TH=10;
c11=[deplonlat(2)-2*TH deplonlat(3)+2*TH];
cmn=[deplonlat(2)+2*TH deplonlat(3)-2*TH];

clf
[ah,ha,H]=krijetem(subnum(3,2));

% Filter ever so slightly before expanding? Don't bother
[r,lon,lat,Plm]=plm2xyz(plmfilt(pr{1},Lup-2,5),degres,[c11 cmn]);
axes(ha(1))
% Plot rotated function on the map
imagefnan(c11,cmn,setnans(r),...
	  kelicol,halverange(r));
[ac,pc(1)]=plotcont(c11+[-2 2],cmn+[2 -2]);
t(1)=title(sprintf('%s = [%i %i %i %i %i %i]/%s','$\mathbf{M}$',...
		   sqrt(2)*M{1},'$\sqrt{2}$'),'interpreter','latex');

for ind=2:length(M)
  r=plm2xyz(plmfilt(pr{ind},Lup-2,5),degres,[c11 cmn],[],[],Plm);
  axes(ha(ind))
  imagefnan(c11,cmn,setnans(r),...
	    kelicol,halverange(r));
  [ac,pc(ind)]=plotcont(c11+[-2 2],cmn+[2 -2]);
  t(ind)=title(sprintf('%s = [%i %i %i %i %i %i]/%s','$\mathbf{M}$',...
		   sqrt(2)*M{ind},'$\sqrt{2}$'),'interpreter','latex');
end

% Cosmetics
nolabels(ah(1:4),1)
nolabels(ha(4:6),2)
set(pc,'linew',0.5)
longticks(ah)
serre(H,1/2,'across')
set(t,'FontS',12)
%moveh(ha(1:3),0.1)
set(ah,'xlim',[c11(1) cmn(1)],'ylim',[cmn(2) c11(2)],...
       'xtick',[80:10:110],'ytick',[-10:10:20])
[bhl,thl]=label(ah,'ll',[],[],[],[],[],[],[]);
deggies(ah)
fig2print(gcf,'tall')
figdisp([],[],[],0)



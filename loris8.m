function loris8(agu)
% loris8(agu)
%
% Synthetic input-output tests with the FISTA inversion
%
% INPUT:
%
% agu     0 Laplacian and FISTA model comparison
%         1 Only FISTA model in one-column format [default]
% 
% Last modified by fjsimons-at-alum.mit.edu, 03/22/2011

defval('agu',1)

cd /u/fjsimons/MyPapers/2011/GJI-2011/GEOsphere

% The one without the holes
fname0='2010-08-19T22h26m01s.mat';

% The one with the holes

% Voxel-based inversion
dirnam='/u/fjsimons/MyPapers/2011/GJI-2011/GEOsphere/pics/';
fname1='2011-02-18T16h32m12s.mat';
voxi=load(fullfile(dirnam,fname1));

% Nope-rather the artifactsy version
%fname1='2011-06-01T12h33m04s.mat';
%fname1='2011-06-01T15h24m18s.mat';
fname1='reconstructions2011-06-01T20h25m39s.mat';
vox2=load(fullfile(dirnam,fname1));

% Wavelet-based inversion
fname2='2011-02-18T17h10m06s.mat';
wavi=load(fullfile(dirnam,fname2));

% Make sure that the operators are identical
difer(strcmp(wavi.operatorfile,voxi.operatorfile)-1,[],[],NaN)
% Make sure that the input files are identical
difer(strcmp(wavi.inputmodelfile,voxi.inputmodelfile)-1,[],[],NaN)

% Map the circles etc
% lat long of centers of holes (lat long)
centers=[30 -90; 0 10; 15 150; -90 0];
% radii of holes (degrees)
radii=[15 20 15 25];
% Make the circles, remember long lat Mollweide
[mx,my]=caploc(centers(:,[2 1]),radii,[],2);
% Make the circles, remember long lat standard
% [lon,lat]=caploc(centers(:,[2 1]),radii,[],1);
% Find the xi,eta coordinates of the circles
% [xi,eta,fid]=sphere2cube(lon,lat);
% Find the indices of the data that fall inside of them
% xic=linspace(-pi/4,pi/4,size(voxi.param,2));
% etac=linspace(-pi/4,pi/4,size(voxi.param,1));
% [xic,etac]=meshgrid(xic,etac);
% This doesn't actually work very well
% for ind=1:6
%   hins{ind}=inpolygon(xic,etac,xi{ind},eta{ind});
% end
[x,y,z]=cube2sphere(7);
% Note that our convention flips later so we prepare this here
for ind=1:6
  x(:,:,ind)=flipud(x(:,:,ind));
  y(:,:,ind)=flipud(y(:,:,ind));
  z(:,:,ind)=flipud(z(:,:,ind));
end
[X,Y,Z]=sph2cart(centers(:,2)*pi/180,centers(:,1)*pi/180,1);
% Shrink the inner radius a bit since IL made his so that the known to be
% zero input really is zero when we calculate the norm inside runplot
il=0.98;
for ind=1:length(radii)
  % Do the dot product and compare to the radius
  hins{ind}=find(180/pi*acos([x(:) y(:) z(:)]*...
                             [X(ind) Y(ind) Z(ind)]')<=il*radii(ind));
end

% Take out the dateline crossings, but not for the polar cap
% Now need to take out the annoying connecting lines
[mxx,myy]=penlift(mx(:,1:3),my(:,1:3));
% Special care for the annoying polar cap
mxx{4}=[mx(:,4) ; -mx(1,4)]; myy{4}=[my(:,4) ; my(1,4)];
mx=mxx; my=myy;

% Make the figure
clf
[ah,ha,H]=krijetem(subnum(3,2));

% Obtain the row sum of the kernels
S=load(voxi.operatorfile);
sumofk=setnans(log10(reshape(full(sum(abs(S.operator),1)),S.boxsize)));

% Color saturation values
dnup=[0.5 2.75];

kelicol

% Plot the sum of the kernels
axes(ah(1))
[cb(1),pp{1},tx{1}]=runplot(ah(1),sumofk,dnup,'log kernel density',...
                      mx,my,hins);
delete(tx{1})

% Obtain the input models
S=load(voxi.inputmodelfile);
inputm=reshape(S.inputparam,S.boxsize);

% Plot the inputmodel
axes(ah(2))
[cb(2),pp{2},tx{2}]=runplot(ah(2),setnans(inputm),[],'input model',...
                      mx,my,hins);

% Obtain output model one
SS=load(voxi.outputname);
outputm=reshape(SS.param,S.boxsize);

% Nope, rather the other one
outputm=vox2.outputmodelCDF42;
outputm=vox2.outputmodelD6withprecond;
outputm=vox2.outputmodelD6noprecond;

% Plot outputmodel one
axes(ah(3))
[cb(3),pp{3},tx{3}]=runplot(ah(3),setnans(outputm),[],'output model 1',...
                      mx,my,hins);

% Obtain output model two
SS=load(wavi.outputname);
outputm=reshape(param2model(SS.param,SS.transformdata),S.boxsize);

% Plot outputmodel two
axes(ah(4))
[cb(4),pp{4},tx{4}]=runplot(ah(4),setnans(outputm),[],'output model 2',...
                      mx,my,hins);

% Plot some of the diagnostic plots, squared misfit normalized by data
% noise versus the l2 norm of the Laplacian of the paramaters
axes(ah(5))
lo(1)=loglog(voxi.l2norm,voxi.chi2,'k+');
hold on
lz(1)=loglog(voxi.l2norm(2),voxi.chi2(2),'o');
le(1)=loglog(voxi.l2norm(end),voxi.chi2(end),'s');
%lh{1}=loglog(voxi.l2norm(100:100:end),voxi.chi2(100:100:end),'+');
hold off
%ylstring1='\textsf{$\|\mathbf{m}-\mathbf{\hat{m}}\|^2_2/\sigma^2/N$}';
ylstring1=sprintf('%s','\chi^2/N misfit norm');
yl(1)=ylabel(ylstring1);
xlstring1='\textsf{$\ell_2$ norm of the model Laplacian}';
xl(1)=xlabel(xlstring1);

% Plot some of the diagnostic plots, squared misfit normalized by data
% noise versus the l1 norm of the paramaters themselves
% Make this the reduced chi2 by dividing by the number of data, not the
% degrees of freedom since we don't pretend we even have this
wavi.chi2=wavi.chi2/length(wavi.noisydata);
axes(ah(6))
lo(2)=loglog(wavi.l1norm,wavi.chi2,'k+');
hold on
lz(2)=loglog(wavi.l1norm(2),wavi.chi2(2),'o');
le(2)=loglog(wavi.l1norm(end),wavi.chi2(end),'s');
%lh{2}=loglog(wavi.l1norm(100:100:end),wavi.chi2(100:100:end),'+');
hold off
ylstring2='\textsf{$\|\mathbf{m}-\mathbf{\hat{m}}\|^2_2/\sigma^2/N$}';
ylstring2=sprintf('%s','\chi^2/N misfit norm');
yl(2)=ylabel(ylstring2);
xlstring2='\textsf{$\ell_1$ norm of the model coefficients}';
xl(2)=xlabel(xlstring2);

% Beautification
shrink(ah(5:6),1.25,1.75)
movev([ah(1:2) cb(1:2)],-.1)
movev([ah(3:4) cb(3:4)],-.01)
movev([ah(5:6)],.07)

set(ah([5 6]),'ylim',[6e3 7e5]/length(wavi.noisydata))
%set(ah([5 6]),'ytick',10.^[4 5])
set(ah([5 6]),'ytick',10.^[0 1],'ygrid','on','yminorgrid','off')
longticks(ah([5:6]))
fs=12;
%set([xl yl],'Interpreter','LaTeX','FontS',fs)
set(xl,'Interpreter','LaTeX','FontS',fs)
set(lz,'MarkerF','w','MarkerE','k','MarkerS',4)
set(le,'MarkerF','w','MarkerE','k','MarkerS',4)
fig2print(gcf,'tall')
set(gcf,'inverthardcopy','off','color','w')

if agu==1
  delete(ha([2 3]))
  delete(cb(3))
  set(ah(1),'position',[getpos(ah(2),1) getpos(ah(1),2:4)])
  set(cb(1),'position',[getpos(cb(2),1) getpos(cb(1),2:4)])
  movev(ha(4:6),-.1)
  movev(cb([2 4]),-.1)
  movev(ah(1),.125)
  movev(cb(1),.125)
  moveh([ah(1) ha(4:6) cb([1 2 4])],-.45)
  set(get(cb(4),'xlabel'),'string','output model')
  movev([ah(1) ha(4:6) cb([1 2 4])],.1)
end

actprint=0;
popts='-zbuffer'',''-r300';
figna=figdisp([],[],popts,actprint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cb,pp,tx]=runplot(ax,data,dnup,cbstring,mx,my,hins)
% sumofk=setnans(sumofk,10);
% Plot IL way
% figure(2); clf ;scalarplotflat(data,'stairs',-0.5); figure(1)
% Plot the circles on top in this space also

% Reorient cubes FJS way
for i=1:6; data(:,:,i)=data(:,:,i)'; end

% Reorder cubes FJS way
data=data(:,:,[5 4 6 2 3 1]);

% Compute the misfit inside of the circles defined by xi and eta
% Check visually that we are looking at the right thing!
for ind=1:4
  incircle=data(hins{ind});
  rmss(ind)=rms(incircle(~isnan(incircle)));
  meen(ind)=mean(incircle(~isnan(incircle)));
  % If they were all NaN's we're in the circle and the norm is 0
  if isnan(rmss(ind)); rmss(ind)=0; end
  if isnan(meen(ind)); meen(ind)=0; end
end
% Now do all but the holes
outcircle=skip(data,cat(1,hins{:}));
rmss(5)=rms(outcircle(~isnan(outcircle)));
meen(5)=mean(outcircle(~isnan(outcircle)));

% Project to Mollweide - this changes the data values
[lonm,latm,data]=cube2moll(data);

% Color saturation values
%defval('dnup',minmax(data(:)));
defval('dnup',halverange(data,65));
% Rounding this symmetric scale to a fifth
dnup=dnup+[-1 1].*mod(dnup(1),.20);
% One-time fix
if dnup==[-0.8 0.8]; dnup=[-0.6 0.6]; end

% Quick look to make sure we're not missing anything here
fprintf(1,'%8.3f\n',rmss)

% mes1=sprintf('rms 1: %5.2f\nrms 2: %5.2f',rmss(1:2));
% mes2=sprintf('rms 3: %5.2f\nrms 4: %5.2f',rmss(3:4));
% mes3=sprintf('rms: %5.2f',rmss(5));
mes1=sprintf('%6.3f%s%4.2f\n%6.3f%s%4.2f',...
             meen(1),'\pm',rmss(1),meen(2),'\pm',rmss(2));
mes2=sprintf('%6.3f%s%4.2f\n%6.3f%s%4.2f',...
             meen(3),'\pm',rmss(3),meen(4),'\pm',rmss(4));
mes3=sprintf('%6.3f%s%4.2f',meen(5),'\pm',rmss(5));
mes4=sprintf('min: %5.2f\nmax: %5.2f',min(data(:)),max(data(:)));

disp(sprintf('Old model range %5.2f to %5.2f',...
             min(data(:)),max(data(:))))
data(data>dnup(2))=dnup(2);
data(data<dnup(1))=dnup(1);
disp(sprintf('New model range %5.2f to %5.2f',...
             min(data(:)),max(data(:))))

% FOR THE PLOTTING ONLY (AFTER THE REPORT) CHANGE THE DATA TO LOOK PRETTIER
data=setnans(data,20);

% Do the actual plotting 
axes(ax)
pc=pcolor(lonm,latm,data); shading flat

% Plot the message
tx(1)=text( 2.25,-1.25,mes1);
tx(2)=text(-2.25,-1.25,mes2);
tx(3)=text( 2.25,1.25,mes3);
tx(4)=text(-2.25,1.25,mes4);
set(tx(1),'horizontala','left')
set(tx([2 4]),'horizontala','right')

% Put on the color bar - not in the function!! messes up
cwidth=0.1875;
cheight=0.01;
pos=[getpos(ax,1)+getpos(ax,3)/2-cwidth/2 ...
     getpos(ax,2)+getpos(ax,4)/10 cwidth cheight];
cb=colorbarf('hor',10,'Helvetica',pos);

set(cb,'xtick',linspace(dnup(1),dnup(2),5),'xlim',dnup)
longticks(cb)
axes(cb)
xlcb=xlabel(cbstring);
axes(ax)

% These are the continents and plate boundaries in mollweide
[a,b]=plotcont([],[],2); 
[c,d]=plotplates([],[],2); 

% These are the circles that are being cut out
hold on
for in=1:length(mx)
  pp(in)=plot(mx{in},my{in},'k-','LineW',1);
end
hold off

% The plotcont following garantees the set(gcf,'NextPlot','add')
axis image off

% Expand the axis slightly, however, might have used XPAND
axis(axis*1.05)


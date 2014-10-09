function loris9(depkm,algo)
% LORIS9(depkm,algo)
%
% This compares the wavelet coefficients by scale and by panel of our
% seismic models as obtained from the breakdown by LORIS5 
%
% d=load('/Users/fjsimons/IFILES/EARTHMODELS/MONTELLI/GNdepths');
% for index=1:length(d); loris9(d(index)/1000,2); end
%
% Last modified by fjsimons-at-alum.mit.edu, 12/02/2011

defval('depkm',474.075);
defval('depkm',2618.700);
defval('depkm',677.250);
defval('depkm',1015.875);
defval('depkm',609.525);
defval('depkm',406.350);
defval('wav','D4')
defval('N',7)
defval('J',N-3-strcmp(wav,'D6'))
defval('fsr',6);
defval('xver',0);
% Regression (1) or not (2)
ops=1;
% Regular least squares (1) or total least squares (2)
defval('algo',2);

% Scale-indexed structure per panel
jscale=cellnan(J+1,129,6);
% Preconditioned it is
precon=[1 1];

% Load the wavelet-transformed pieces of the model
% Montelli
[vw1,depkm]=loris5('D4',7,4,precon,depkm,0,1);
% Ritsema
[vw2,depkm]=loris5('D4',7,4,precon,depkm,0,2);

% Identify the scales in the mr=1 wavelet transform
% Do one more than J so the last one is the scaling coefficients 
vwlev=cube2scale(N,[J J]+1,1);

% This'd be the legend for the chunk numbers ordered sequentially
% Remember the reordering by the others
pv=[5 4 6 2 3 1];
legs={'1 Pacific','2 Antarctica','3 Asia',...
      '4 S America','5 N America','6 Africa'};
legs={legs{pv}};

% Are we still talking about the same thing? 
% Quick plot of the wavelet coefficients
if xver==1
  % Now calculate the inverse wavelet transforms
  switch wav
   case 'D2'
    % The 2-tap Haar/Daubechies wavelet transform
    vwrec1=angularD2WT(vw1,[J J],'inverse',1);
    vwrec2=angularD2WT(vw2,[J J],'inverse',1);
   case 'D4'
    % The 4-tap Daubechies wavelet transform
    vwrec1=angularD4WT(vw1,[J J],precon,'inverse',1);
    vwrec2=angularD4WT(vw2,[J J],precon,'inverse',1);
   case 'D6'
    % The 6-tap Daubechies wavelet transform
    vwrec1=angularD6WT(vw1,[J J],precon,'inverse',1);
    vwrec2=angularD6WT(vw1,[J J],precon,'inverse',1);
  end
  % Panel legend verification
  for pans=1:6
    vwpan=vwrec1; vwpan(:,:,skip(1:6,pans))=0;
    clf
    [cb,xcb]=quickplot(N,J,vwpan,depkm);
    title(legs{pans})
    pause
  end

  % Make a plot of either model or their wavelet coefficients
  clf
  [cb,xcb]=quickplot(N,J,vwrec1,depkm);
  pause
  clf
  [cb,xcb]=quickplot(N,J,vwrec2,depkm);
  pause
  clf
  [cb,xcb]=quickplot(N,J,vw1,depkm);
  pause
  clf
  [cb,xcb]=quickplot(N,J,vw2,depkm);
end

% Start figure
clf
[ah,ha,H]=krijetem(subnum(J+1,6));
ja=zeros(prod(size(ah),1));
axl=[-1 1 -1 1]*1.1;
% Line style for fit and its uncertainty
colx='k-'; colxx='k--';
Rav=0;

% Make a scatter plot of the coefficients by scale
for jndex=1:J+1
  logs{jndex}=sprintf('scale %i (wavs)',jndex);
  for nchunk=1:6
    axes(ah(nchunk+(jndex-1)*6))
    vw1J=vw1(:,:,nchunk);
    vw2J=vw2(:,:,nchunk);
    % Identify where the scales are located per chunk
    vw1J=vw1J(vwlev==jndex);
    vw2J=vw2J(vwlev==jndex);

    % Make it all relative to the overall max
    maxx=max(abs([vw1J ; vw2J]));
    % Use same benchmark or else lose slope
    vw1J=vw1J/maxx;
    vw2J=vw2J/maxx;

    % Make the 2D histograms
    [Hi,c11,cmn,HH,ybine]=bindens(vw1J,vw2J,15,15);
    
    switch ops
     case 1
      % Do the regression
      cmap=gray;
      plotthis=decibel(100*HH/sum(HH(:)),100);
      plotthis(isinf(plotthis))=NaN;
      % Cut this off below the one percent level
      % But make the maximum value really black
      cax=[-20 max(plotthis(:))];
      imagefnan(c11,cmn,plotthis,cmap,cax,[],1,0)

      switch algo
       case 1
	% Calculate the regression line
	[P,S]=polyfit(vw1J,vw2J,1);
	CP=(inv(S.R)*inv(S.R)')*S.normr^2/S.df;
	% Prepare the report
	tx2=sprintf('S/P %4.2f%s%4.2f',P(1),'\pm',sqrt(CP(1)));
	Pp=P; Pm=P;
	% Bit of a half-assed measure but all right
	Pp(1)=P(1)+sqrt(CP(1));
	Pm(1)=P(1)-sqrt(CP(1));
	% Don't use the "prediction" error, no fun
	PX=polyval(P,axl(1:2));
	% But predict based on the "functional" error
	PXp=polyval(Pp,axl(1:2));
	PXm=polyval(Pm,axl(1:2));
       case 2
	% Calculate the regression line
	P=tls(vw1J,vw2J);
	PX=polyval(P,axl(1:2));
	% Prepare the report
	tx2=sprintf('S/P %4.2f',P(1));
      end
      % Calculate the correlation coefficient
      [R,RP]=corrcoef(vw1J,vw2J);
      % Plot regression if correlation is significant
      % Plot regression if correlation is big
      % Clean up
      axis equal tight 
      axis(axl)
      longticks(gca,1/2)
      set(gca,'xgrid','off','ygrid','off')
      hold on
      plot(axl(1:2),[0 0],':')
      plot([0 0],axl(3:4),':')
      % If high and significant
      if R(2)>=0.35 && RP(2)<0.05
	p1=plot(axl(1:2),PX,colx);
	try
	  p2=plot(axl(1:2),PXp,colxx);
	  p3=plot(axl(1:2),PXm,colxx);
	end
	% On laptop
	%[b2,t2]=boxtex('ur',gca,tx2,fsr,[],[],0.7+0.2*[algo==2]);
	% On desktop
	[b2,t2]=boxtex('ur',gca,tx2,fsr,[],0.5,0.5+0.2*[algo==2]);
	% This doesn't work so well
	top(b2,gca)
	% Keep track of which ones have had this done and which not
	ja(nchunk+(jndex-1)*6)=1;
      end
      % If significant
      if RP(2)<0.05
	tx1=sprintf('R %4.2f',R(2));
	% On laptop
	%[b1,t1]=boxtex('ll',gca,tx1,fsr,[],[],0.9);
	% On desktop
	[b1,t1]=boxtex('ll',gca,tx1,fsr,[],0.5,0.7);
	% This doesn't work so well
	top(b1,gca); 
	% Keep track of which ones have had this done and which not
	ja(nchunk+(jndex-1)*6)=1;
      end
      hold off
     case 2
      % Simply plot the data
      plot(vw1J,vw2J,'y.','MarkerS',2)
    end
    % Titles and labels
    if jndex==1
      ts(nchunk)=title(legs{nchunk});
    end
    if nchunk==1
      tc(jndex)=ylabel(logs{jndex});
    end
  end
  % Keep the average of the absolute value
  Rav=Rav+abs(R(2));
end

% Make the average
Rav=Rav/(J+1)/6;

disp(sprintf('Average abs(R) is %4.4f at %i km depth',Rav,round(depkm)))

set(findobj(tc,'string',sprintf('scale %i (wavs)',J+1)),...
    'string',sprintf('scale %i (scals)',J)) 

% Various cosmetics
nolabels(ha(J+2:end),2)
nolabels(ah(1:end-6),1)

axes(H(J+1,3))
yll(1)=xlabel(sprintf(...
    'Montelli (2006) %s wave model at %i km',...
    '{\itP}',round(depkm)));
moveh(yll(1),1.5)

axes(H(J/2+1,6))
yll(2)=ylabel(sprintf(...
    'Ritsema (2010) %s wave model at %i km',...
    '{\itS}',round(depkm)));
set(gca,'YaxisLoc','r')
set(yll(2),'Rotation',-90)
moveh(yll(2),.25)

%set(ah,'CameraV',6.5)
%set(ah(logical(ja)),'CameraV',6.5)
%set(ah(~logical(ja)),'CameraV',5.75)
% For SPIE
set(ah(logical(ja)),'CameraV',4.4)
set(ah(~logical(ja)),'CameraV',4.0)
movev(yll(1),.5)
% End for SPIE

movev(ah,0.01)
fig2print(gcf,'portrait')
actprint=0;
figna=figdisp([],sprintf('%4.4i_%i',round(depkm),algo),[],actprint);
system(sprintf('epstopdf %s.eps',figna));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cb,xcb,pgw]=quickplot(N,J,vw,depkm,dax)
% Font size of the label
defval('fs',8)
% Color map and saturation as percentiles
defval('colmap','kelicol');
% Grid information
wg.N=N; wg.J=J;
% Not at ALL the same as using HALVERANGE
colperc=[1 99];
% May go with HALVERANGE after all
defval('dax',round(halverange(vw,50,NaN)));
% Explicit and absolute color limits of the VALUE of the coeffs
defval('dax',prctile(vw(:),colperc));

% The actual plotting
% Cancel the grid
wg=[];
[~,~,~,pgw]=plotoncube(vw,'2D',1,[],[],[],dax,[],0,100,wg);
plotcont([],[],9)

% Color bar etc
colpos=[0.5616    0.1714+0.025    0.3143    0.0298];
[cb,xcb]=addcb(colpos,dax,dax,colmap,range(dax)/4);
set(cb,'fonts',fs)
set(xcb,'string',sprintf(...
    'wavelet coefficients at %i km',round(depkm)),'fonts',fs)

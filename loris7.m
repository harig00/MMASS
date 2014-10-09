function loris7
% LORIS7
%
% This is showing the breakdown by scale and by panel of our seismic
% models as obtained from the single-chunk breakdown by LORIS5 
%
% Last modified by fjsimons-at-alum.mit.edu, 02/17/2011

defval('wav','D4')
defval('N',7)
defval('J',N-3-strcmp(wav,'D6'))
defval('ndeps',129)
defval('fs',11);

% Scale-indexed structure per panel
jscale=cellnan(J+1,129,6);

for panel=1:6
  % Load the data
  [~,~,vwstats1]=loris5('D4',7,4,[1 1],[],0,1,panel);
  [~,~,vwstats2]=loris5('D4',7,4,[1 1],[],0,2,panel);
  % Extract by scale
  for scal=1:J+1
    jscale1{scal}(:,panel)=cat(1,vwstats1.normsq{:,scal});
    jscale2{scal}(:,panel)=cat(1,vwstats2.normsq{:,scal});
    % Keep track for the axes
    maxi(scal)=max([jscale1{scal}(:) ; jscale2{scal}(:)]);
  end
end

% Make sure that the sum of the proportions per chunk sums to 100
disp(sprintf('Norm squared percentages add up to about %6.2f %s',...
             mean(sum([jscale1{:}],2)),'%'))

% Get the common depth axis
deps=cat(1,vwstats1.deptkm{:});
deps=[vwstats1.deptkm{:}];

% Prepare stuff common to all axes
ytix=[0 410 660 1000 2000 2899];
ytil={'0' '410' '660' '1000' '2000' 'CMB'};

% Start figure
clf
[ah,ha,H]=krijetem(subnum(2,J+1));
xmin=[0.4 1.25 4 8.1 40];
xmin=zeros(J+1,1);

% Now plot this business
for scal=1:J+1
  xmax=max(ceil(maxi(scal)*10)/10,xmin(scal));
  % For the Montelli model
  axes(ah(scal))
  p1{scal}=plot(jscale1{scal},deps);
  hold on
  xlim([0 xmax])
  g1{scal}=plot(xlim,[ytix(2:end-1)' ytix(2:end-1)'],'k:');
  hold off
  t1(scal)=title(sprintf('scale %i (wavs)',scal));
  % y1(scal)=ylabel('depth (km)');
  % For the Ritsema model
  axes(ah(scal+J+1))
  % hold on
  p2{scal}=plot(jscale2{scal},deps);
  hold on
  xlim([0 xmax])
  g2{scal}=plot(xlim,[ytix(2:end-1)' ytix(2:end-1)'],'k:');
  hold off
  % y2(scal)=ylabel('depth (km)');
end

% Cosmetics
cols={'k' 'y' 'r' 'b' 'm' 'g'};
lins={'-' '-' '--' '--' '-' '--'};
for ind=1:6
  set(rindeks([p1{:}],ind),'Color',cols{ind},'lines',lins{ind})
  set(rindeks([p2{:}],ind),'Color',cols{ind},'lines',lins{ind})
end

% Remember the reordering by the others which is reported in the paper
pv=[5 4 6 2 3 1];
axes(ah(J+2))
lel=legend({'5','4','6','2',...
          '3','1'},'Location','SouthEast');

set(findobj(t1,'string',sprintf('scale %i (wavs)',J+1)),...
    'string',sprintf('scale %i (scals)',J)) 
set([p1{:}],'linew',1)
set([p2{:}],'linew',1)
set(ah,'ydir','rev')
set(ah,'ytick',ytix,'ytickl',ytil,'ylim',minmax(ytix))
%delete(y1(2:end))
%delete(y2(2:end))
longticks(ah)
nolabels(ha(3:end),2)
noticks(ah(1:J+1),1)
shrink(ah,.8,1)
fig2print(gcf,'portrait')

axes(ah(J+1+ceil((J+1)/2)))
xlstring2=['\textsf{contribution to $\ell_2$ model norm per scale (\%)}'];
xl(1)=xlabel(xlstring2);
set(xl(1),'Interpreter','LaTeX','FontS',fs+1)

serre(H',1/2,'down')

axes(ah(J+1))
yll(1)=ylabel('Montelli (2006) {\itP} wave model');
axes(ah(end))
yll(2)=ylabel('Ritsema (2010) {\itS} wave model');
set(ah([J+1 end]),'yaxisl','r');
set(yll,'rotation',-90)
moveh(yll,2)
figdisp()


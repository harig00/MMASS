function swvals2d
% SWVALS2D
% 
% FIGURE 3 of SIMONS & WANG
%
% Makes a plot of the eigenvalues of the concentration problem to the
% Cartesian disk while differentiating the orders.
%
% SEE ALSO:
%
% SWFRIED2D, SDWVALS
%
% Last modified by fjsimons-at-alum.mit.edu, 04/24/2009

% Prepare for plotting
clf
[ah,ha,H]=krijetem(subnum(2,2));
yls=[-0.1 1.1];
symbs={'o','x','s','+','v','*','^','d','<','>','p','h',...
       'o','x','s','+','v','*','^','d'};
xmax=60;
NN=[3 11 24 42];

% Do this for a variety of Shannon numbers
for index=1:length(NN)
  clear E EV EM

  N=NN(index);
  
  % The upper of number of orders should "cover" all the cases down to
  % very low values
  M=11;
  
  % Calculate the radial functions but only care about the eigenvalues
  for m=0:M
    [E,EV{m+1}]=swdisk(m,N,2*N,[],NaN,'DV');
    EM{m+1}=repmat(m,1,2*N);
  end
  % Sort them all according to their eigenvalue
  EV=[EV{:}]; EM=[EM{:}];
  [EV,i]=sort(EV,'descend'); EM=EM(i);
  % Repeat the nonzero orders twice and take as many as you had
  dbl=~~EM+1;
  % The indexing sequence in the non-repeated vectors
  seq=gamini(1:length(EV),dbl);
  % Make the eigenvalue sequence with the repeats
  EV=gamini(EV,dbl);
  % Make the order sequence with the repeats
  EM=EM(seq);
  
  % Now do the plots
  axes(ah(index))
  for ondi=1:length(EV)
    p(ondi,index)=plot(ondi,EV(ondi),symbs{EM(ondi)+1});
    hold on
  end
  hold on
  plot([N N],yls,'k:')
  plot([0 xmax],[0.5 0.5],'k:')
  plot([0 xmax],[0 0],'k:')
  plot([0 xmax],[1 1],'k:')
  set(ah(index),'xlim',[0 xmax],'ylim',yls,'xgrid','off','ygrid','off',...
		'xtickl',[1 10:10:xmax],'xtick',[1 10:10:xmax],...
		'ytick',[0:0.25:1])
  drawnow
end

% Save the cosmetics for the very end
% Now make the plot beautiful
axes(ah(4))
fb=fillbox([2 18 0.88 -0.05],'w');
for ondi=1:12
  ypo=0+0.075*(ondi-1);
  pl(ondi,1)=plot(4,ypo,symbs{ondi});
  hold on
  tl(ondi,1)=text(7,ypo,sprintf('m = %s %i','\pm',ondi-1),'FontS',8);
end

longticks(ah)
set([p(~~p(:)) ; pl(~~pl(:))],'MarkerS',4,'MarkerF',grey,'MarkerE','k')
axes(ha(1))
al(1)=ylabel('eigenvalue \lambda');
axes(ha(2))
al(2)=ylabel('eigenvalue \lambda');
axes(ha(2))
xl(1)=xlabel('rank');
axes(ha(4))
xl(2)=xlabel('rank');

nolabels(ha(3:4),2)
nolabels(ah(1:2),1)

serre(H',1/2,'down')
serre(H,1/2,'across')

for ind=1:4
  xx(ind)=xtraxis(ah(ind),NN(ind),sprintf('N2D = %i',NN(ind)));
end
longticks(xx)

set([xl al],'FontS',13)
set([ ah],'FontS',12)

fig2print(gcf,'portrait')

figdisp([],[],[],1)

function swfried2d(method)
% swfried2d(method)
%
% FIGURE 4 of SIMONS & WANG, GEM 2011
%
% INPUT:
%
% method   'SE' by Slepian "extension" [faster!]
%          'GL' by direct Gauss-Legendre integration [slow!]
%          'RS' by direct Riemann summation [slow!]
%          'compare' by comparison of the methods
%
% Makes a fried-egg plot for +- angular orders for concentration to a
% disk-shaped domain. If 'method' is 'GL' The eigenvalues will be ever so
% slightly different from the "truth" that we know so well - they may not
% occur exactly in pairs, etc. But overall this is pretty well behaved to
% a part in a million or so starting from reasonable choice.
%
% SEE ALSO:
%
% SWVALS2D
%
% Last modified by fjsimons-at-alum.mit.edu, 03/03/2011

% eps2png -width 850 swfried2d.eps 
% convert swfried2d.png swfried2d.gif

% Somehow this performs best with a KEYBOARD before the print statement

% Here we both simultaneously check and compute the eigenfunctions

defval('method','SE');

% Make a circle of some description
circn=256;
XY=[cos(linspace(0,2*pi,circn)) ; sin(linspace(0,2*pi,circn))]';

% What is the mean distance between the points
mdis=mean(sum(sqrt(diff(XY).^2),2));

% Perform the localization
N=42; % The Shannon number 
J=30; % The number of Slepian functions that will be calculated

% Define the (resolution of the) output grid
j=7;
xp=linspace(-2,2,2^j);
yp=linspace(-2,2,2^j);
[XP,YP]=meshgrid(xp,yp);
XYP=[XP(:) YP(:)];

% Now do the actual calculation
switch method
  case {'GL','RS'}
   % This returns all of the numerical eigenvalus
   [G,H,V,K,XYP,XY,A]=localization2D(XY,N,J,method,[],[],[],XYP);
   G=reshape(G,2^j,2^j,J);
 
   % Make the figure - with G/E or GV/EV for comparison
   theG=G; 
   theV=V;
 case 'SE'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we need to compare this to the other way of calculating it
  % The upper of number of orders should "cover" all the cases above
  M=8;
  % Calculate at least all of those that have eigenvalues of larger than 0.5
  [Nm,Nsum]=nsubm(N,M);
  % Form the spatial eigenfunction
  R=sqrt(XP.^2+YP.^2);
  TH=atan2(YP,XP);
  [Rx,iR,jR]=unique(R);
  % Calculate the radial functions at the unique radial points that will be needed
  for m=0:M
    [E{m+1},EV{m+1}]=swdisk(m,N,ceil(Nm(m+1)),[],Rx,'SE',0);
    EM{m+1}=repmat(m,1,ceil(Nm(m+1)));
  end
  % Sort them all and compare eigenfunctions with the same eigenvalue
  E=[E{:}]; EV=[EV{:}]; EM=[EM{:}];
  [EV,i]=sort(EV,'descend'); E=E(:,i); EM=EM(i);
  % Repeat the nonzero orders twice and take as many as you had
  dbl=~~EM+1;
  % The indexing sequence in the non-repeated vectors
  seq=indeks(gamini(1:length(EV),dbl),1:J);
  % This will change the sine into the cosine later
  % If there is no difference, it is a pair
  addph=~diff([0 seq]); 
  % Make the eigenvalue sequence with the repeats
  EV=indeks(gamini(EV,dbl),1:J);
  % Reassemble and combine with the sine/cosine functions
  GV=nan(size(R,1),size(R,2),J);
  for i=1:length(EV)
    % Just leave it here - it's the GL that is in trouble
    GV(:,:,i)=reshape(E(jR,seq(i)),size(R)).*...
	      cos(EM(seq(i)).*TH-addph(i)*pi/2);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Make the figure - with G/E or GV/EV for comparison
  theG=GV;
  theV=EV;
 case 'compare'
  keyboard
  % Not properly written but largely irrelevant now
  % Check how close these eigenvalues are now
  difer(EV-V,4)
  % Also compare the eigenfunctions and make sure that the normalization
  % is identical in both cases.
end
  
% They should differ slightly - in the sign - in the phase
clf
[ah,ha,H]=krijetem(subnum(6,5));

c11=[min(XYP(:,1)) max(XYP(:,2))];
cmn=[max(XYP(:,1)) min(XYP(:,1))];

defval('bw',1)

for index=1:length(ah)
  axes(ah(index))
  % Check if the color saturation is appropriate
  if bw==1
    imagefnan(c11,cmn,setnans(theG(:,:,index)),gray(21),[-1 1],grey(5))
  else
    imagefnan(c11,cmn,setnans(theG(:,:,index)),kelicol,[-1 1])
  end
  hold on
  % Plot the circle
  plot(XY(:,1),XY(:,2),'k--')
  % Plot the square
  plot(c11([1 2 2 1 1]),cmn([1 1 2 2 1]),'k-')
  switch method
   case {'GL','RS'}
    tl(index)=title(sprintf('%s_{%i} = %9.6f','\lambda',...
			    index,theV(index)));
   case 'SE'
    tl(index)=title(sprintf('%s_{%i} = %9.6f ; m = %i',...
			    '\lambda',...
			    index,theV(index),...
			    (-1).^[addph(index)-1]*EM(seq(index))));
  end
  axis off
  drawnow
end

% Cosmetics (best after KEYBOARD and manually line by line)
movev(tl,-0.3)
set(tl,'FontS',8)
serre(H',[],'down')

fig2print(gcf,'tall')
figdisp([],[],[],0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=setnans(data)
data(abs(data)<max(abs(data(:)))/100)=NaN;


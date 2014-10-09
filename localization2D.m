function varargout=localization2D(XY,N,J,method,Nqx,Nqy,dxdy,XYP)
% [G,H,V,K,XYP,XY,A,c11,cmn]=LOCALIZATION2D(XY,N,J,method,Nqx,Nqy,dxdy,XYP)
% 
% Calculates spatiospectrally concentrated Slepian functions in the
% two-dimensional Cartesian plane at the desired spatial points, for an
% abitrarily shaped region. These are windows to be applied to certain
% data of which you know the specific coordinates or the x and y
% spacings. If you don't, put in a reasonable sampling interval in km.
%
% INPUT:
%
% XY       [X(:) Y(:)] coordinates of a closed curve [default: circle]
% N        Shannon number [default: 20]
% J        Number of eigentapers [default: 2*N]
% method   'GL' by Gauss-Legendre integration [default]
%          'RS' by equal-area Riemann summation
% Nqx      Number of x quadrature points [dx-based for 'RS' / 32 for 'GL']
% Nqy      Number of y quadrature points [dy-based for 'RS' / 32 for 'GL']
% dxdy     [dx dy] the x and y spacing for the return grid OR
% XYP      [X(:) Y(:)] evaluation points [default: rectangle inscribing region]
% c11,cmn  Coordinates of the grid corner points for plotting
%
% OUTPUT:
%
% G        The first J space-concentrated eigentapers at XYP [size(XYP) J]
% H        The first J spacelimited eigentapers at XYP [size(XYP) J]
% V        The eigenvalues of the tapers [vector]
% K        The circular bandlimit wavenumber for this problem
% XYP      Return coordinates for H and G [unwrapped two-column matrix]
% XY       The spatial curve - again
% A        The approximate area of the domain enclosed by the points
%
% EXAMPLES:
%
% localization2D('demo1'); % For the circle
% localization2D('demo2'); % For a random blob
% localization2D('demo3'); % For a profile through the circle
% localization2D('demo4'); % For a profile through the blob
% localization2D('demo1',[],[],'RS'); 
% localization2D('demo1',[],[],'GL'); 
% localization2D('demo2',[],[],'RS'); 
% localization2D('demo2',[],[],'GL'); 
% localization2D('demo4',[],[],'RS');
% localization2D('demo4',[],[],'GL');
% localization2D('demo5') % Compare eigenvalues for circle
% localization2D('demo6') % Compare eigenvalues for blob
% localization2D('demo7') % For a square
%
% SEE ALSO:
%
% SWDISK, KERNELC2D, LOCALIZATION, KERNELC, DEVILLIERS
%
% Last modified by dongwang-at-princeton.edu, 05/14/2010
% Last modified by fjsimons-at-alum.mit.edu, 06/17/2010

% Supply defaults - make sure to hit the diameter exactly
defval('circn',41)
% Note that the discretization of the boundary curve makes for the
% integration accuracy of the circle in the Gauss-Legendre case, and not
% so much the number of nodes
defval('XY',[cos(linspace(0,2*pi,circn)) ; sin(linspace(0,2*pi,circn))]')

if ~isstr(XY)
  defval('N',20)
  defval('method','GL');

  % Make sure the XY of the curve has two columns
  if ~isempty(find(size(XY,2)==2))
    if size(XY,1)<size(XY,2)
      XY=XY';
    end
  else
    error('Coordinates of the curve not as expected')
  end

  % Find limits in x and y
  xlimt=minmax(XY(:,1)); xrange=xlimt(2)-xlimt(1);
  ylimt=minmax(XY(:,2)); yrange=ylimt(2)-ylimt(1);

  % Should be roughly compatible with the resolution of the curve but
  % let's say we're going to assume units of km coming in - thus 5 km.
  % Got to change this if you're not happy! But too low and you might be waiting a long while 
  defval('dxdy',[5 5]/100)
  disp(sprintf('Resolution %gx%g',dxdy(1),dxdy(2)))

  % Find the integration nodes and the corresponding weights
  switch method
   case 'RS'  % The equal-area Riemann summation method
     % Construct a grid with the region inscribed
     % How fine will this grid be? The default will be as fine as the data,
     % thus the quadrature grid equals the data grid
     defval('Nqx',xrange/dxdy(1)) % No rounding, purely formal here
     defval('Nqy',yrange/dxdy(2)) % No rounding, purely formal here
     % So if Nqx and Nqy are blank, quadrature-dx and q-dy are data-dxdy
     qdx=xrange/Nqx;
     qdy=yrange/Nqy;
     % Now make the actual quadrature-grid points: the full mesh set of
     % possible node points in this quadrature grid (q-grid). Watch the
     % upper limit of this range, go all the way up to the boundary
     qx=xlimt(1)+qdx/2:qdx:xlimt(2); 
     qy=ylimt(1)+qdy/2:qdy:ylimt(2); 
     % Report on the method
     disp(sprintf('\nRiemann summation method on a %i x %i grid\n',...
		  length(qy),length(qx)))
     % This is the quadrature grid
     [QX,QY]=meshgrid(qx,qy);
     % The midpoint indices of this subset that fall inside of the region...
     QinR=find(inpolygon(QX,QY,XY(:,1),XY(:,2)));
     % ...and the corresponding weights; this is a column vector
     w=repmat(qdx*qdy,length(QinR),1);
   case 'GL' % The Gauss-Legendre method
    % Default number of GL node points
    defval('Nqx',32);
    defval('Nqy',32);
    % Report on the method
    disp(sprintf('\nGauss-Legendre method with %i x %i integration nodes\n',...
	 Nqy,Nqx))
    % Calculate "blanks", x-axis GL points regardless of the data range
    [wqx,qx]=gausslegendrecof(2*Nqx-1);

    % Calculate y-axis GL points spanning the entire actual y-range;
    [wqy,qy]=gausslegendrecof(2*Nqy-1,[],ylimt);
    
    % For each of these points in y, find the appropriate ranges of x's
    [xints,yp,xp,forreal]=phicurve([XY(:,2) XY(:,1)],qy);
    
    % Check out you're hatching well by doing this:
    % plot(XY(:,1),XY(:,2),'b'); hold on; plot(xp,yp,'k'); axis image

    % For each of these intervals, calculate the required scalings
    % How many intervals total are we talking about?
    % Nall=sum(sum(~~xints(:,1:2:end),1));
    Nall=sum(forreal(:))/2;
    % Initialize the quadrature points
    [QX,QY,w]=deal(nan(Nqx,Nall));
    
    % From GAUSSLEGENDRE: this rescales the nodes to the correct interval
    % Could probably cut this loop, too, by making them all at once,
    % leave for later.
    Nrun=0;
    for yindex=1:length(qy)
      % How many horizontal intervals are we talking about at this y?
      % Nxint=sum(~~xints(yindex,:))/2;
      Nxint=sum(forreal(yindex,:))/2;
      abset=reshape(xints(yindex,1:Nxint*2),2,Nxint);
      % Beginning and end points of those intervals
      a=abset(1,:);
      b=abset(2,:);
      % Now produce the scaled nodes, a column per interval
      exs=repmat(a,length(qx),1)+(qx+1)/2*(b-a);
      % And produce the weights that go with it
      % And multiply this with the current y-node and weight
      % And produce copies of the y-nodes, as many as needed
      % And put them into two big vectors of weights and points
      QX(:,Nrun+1:Nrun+Nxint)=exs;
      QY(:,Nrun+1:Nrun+Nxint)=repmat(qy(yindex),size(exs));
      % Note that w is a matrix here
      w(:,Nrun+1:Nrun+Nxint)=wqx*(b-a)/2*wqy(yindex);
      Nrun=Nrun+Nxint;
    end
    
    % All of these are automatically good, but must be a column
    QinR=[1:length(QX(:))]';
  end

  % Done with the discretization of the integral kernel. Only now do you
  % know how big this will be, so now provide default on the number of
  % eigenvalues and protect against having perhaps wanted too many
  defval('J',2*N)
  if J>length(QinR)
    error('Increase quadrature resolution if you want this many eigenvalues')
  end

  % Report on the area using these approximations
  A=sum(w(:));
  disp(sprintf('Approximate area is %8.3f',A))
  
  % The wavenumber of the circular bandlimitation
  K=sqrt(4*pi*N/A);

  disp(sprintf('Starting kernel calculation of size %i x %i',...
	       length(QinR),length(QinR)))
  % Now calculate the integral kernel at the same unwrapped set of
  % quadrature points in both of its slots
  D=kernelc2D([QX(QinR) QY(QinR)],[QX(QinR) QY(QinR)],K);
  
  % Square root of weight matrix - does precomputing it actually save time?
  sqrtW=sqrt(diag(w(:)));

  % Calculate the eigenvectors, [i.e. G and H on the q-grid] and eigenvalues
  % Here I'd substitute EIGS in case you don't want the whole set
  t0=clock;
  if J<length(D)
    OPTS.disp=0;
    [QGH,V]=eigs(sqrtW*D*sqrtW,J,'LM',OPTS);
  else
    [QGH,V]=eig(sqrtW*D*sqrtW);
  end
  disp(sprintf('while EIG[S] took %8.4f s',etime(clock,t0)))

  % Sort eigenvalues ; protect against small complex values
  [V,isrt]=sort(sum(real(V)),'descend');

  % Sort and rescale eigenvectors THIS SEEMS NOT KOSHER NEED TO CHECK
  % Don't use inv, as it is slow - same as in SDWCAPT
  QGH=QGH(:,isrt)./repmat(diag(sqrtW),1,size(QGH,2));

  % Restrict to the number of wanted ones - note that the number you had
  % calculated depended on the size of the quadrature grid
  QGH=QGH(:,1:J); V=V(1:J);
  
  % Test orthonormality of this subset; one check for G and H as we haven't
  % rescaled anything yet - this is the integral over the quadrature grid
  % So given that the quadrature is over the region only, this only checks
  % the potential orthonormalization of the H's; to check the G's we'd have
  % to perform the quadrature over the entire region, which can never be big
  % enough to be sure. Only work with the chosen eigenvalues - if there
  % is a large null-space, you'll pick up lots of inaccuracy here
  orv1=QGH'*diag(w(:))*QGH; orv2=orv1-diag(diag(orv1));
  difer(diag(orv1)-1,[],[],'Orthonormality  (on-diagonal) checks out')
  difer(orv2,[],[],'Orthonormality (off-diagonal) checks out')
  % If it doesn't, sometimes, it could be due to numerical degeneracy on
  % the eigenvalues. Let's just live with it. Or to the small size of the
  % calculation domain in comparison with the concentration domain.
  
  % Make sure that this normalization is consistent with SWDISK so that
  % we can compare later on.

  % Now where do you want these eigenfunctions really evaluated? Either at
  % a specific set of points, or else simply at the data grid.
  defval('XYP',NaN)
  if isnan(XYP)
    special=0;
    % Now you're not wanted anywhere special, thus you get the function
    % evaluated on the region-inscribing rectangular evaluation-grid with
    % the data spacing. But save yourself the trouble of redoing this if
    % you had accepted all the defaults so far and done Riemann
    if ~strcmp(method,'RS') || Nqx~=xrange/dxdy(1) || Nqy~=yrange/dxdy(2)
      % The evaluation grid will be the data-grid determined by dxdy
      %ex=xlimt(1)+dxdy(1)/2:dxdy(1):xlimt(2);
      %ey=ylimt(1)+dxdy(2)/2:dxdy(2):ylimt(2);
      % Actually, try to make the default a bigger region now
      ex=xlimt(1)+dxdy(1)/2-xrange/2:dxdy(1):xlimt(2)+xrange/2;
      ey=ylimt(1)+dxdy(2)/2-yrange/2:dxdy(2):ylimt(2)+yrange/2;
      % The evaluation grid (e-grid)
      [EX,EY]=meshgrid(ex,ey);
    else
      % The evaluation grid (e-grid) is the quadrature grid (q-grid)
      EX=QX; EY=QY;
    end
    % Size of rectangular matrix the same size as the
    % meshgridded evaluation-grid (based on the quadrature grid or not)
    GHshape=[size(EX) J];
    % The output evaluation grid, unwrapped, as needed for KERNELC2D
    XYP=[EX(:) EY(:)];
    if ~strcmp(method,'RS') || Nqx~=xrange/dxdy(1) || Nqy~=xrange/dxdy(2)
      % These now are the requested points in this case
      % The indices of this new subset that fall inside of the region
      EinR=find(inpolygon(EX,EY,XY(:,1),XY(:,2)));
    else
      EinR=QinR;
    end
  else
    special=1;
    % Work with the points that you had requested, make sure they are also
    % unwrapped and thus have two columns
    if ~isempty(find(size(XYP,2)==2))
      if size(XYP,1)<size(XYP,2)
	XYP=XYP';
      end
    else
      error('Output coordinates not as expected')
    end
    % Find the subset of requested points that's inside of the region
    EinR=find(inpolygon(XYP(:,1),XYP(:,2),XY(:,1),XY(:,2)));
    % Size of matrix with two columns te same size as the
    % requested matrix of coordinates
    GHshape=[size(XYP) J];
  end

  % Calculate the values at the wanted points
  disp(sprintf('Starting kernel calculation of size %i x %i',...
	       size(XYP,1),length(QinR)))
  Di=kernelc2D([XYP(:,1) XYP(:,2)],[QX(QinR) QY(QinR)],K);

  % This then, evaluated on the requested (list or grid of) points, is G 
  G=Di*diag(w(:))*QGH*diag(1./V);
  % And this is the equivalent H, which is zero outside of the region
  H=zeros(size(G));
  H(EinR,:)=G(EinR,:);

  % Final calculation
  if special==0
    % G and H have the size of the evaluation-grid
    G=reshape(G,GHshape);
    H=reshape(H,GHshape);
  end

  % Define the plotting grid for use in e.g. IMAGEF
  c11=[min(XYP(:,1)) max(XYP(:,2))];
  cmn=[max(XYP(:,1)) min(XYP(:,2))];
  
  % Make output - or suppress it
  varns={G,H,V,K,XYP,XY,A,c11,cmn};
  varargout=varns(1:nargout);
elseif strcmp(XY,'demo1')
  defval('method','RS')
 [G,H,V,K,XYP,XY]=localization2D([],[],[],method);
 demoplots1(G,H,V,K,XYP,XY)
elseif strcmp(XY,'demo2')
  [x,y]=blob;
  [x,y]=randcirc([],[],[],1);
  XY=[x(:) y(:)];
  defval('method','RS')
  [G,H,V,K,XYP,XY]=localization2D(XY,[],[],method);
  demoplots1(G,H,V,K,XYP,XY)
elseif strcmp(XY,'demo3')
  defval('method','RS')
  % Get the functions on the inscribed data grid
  [G1,H1,V1,K1,XYP1,XY]=localization2D([],[],method);
  % Now find somewhere specific
  wheir=round(size(XYP1,1)/2);
  wheir=ceil(rand*size(XYP1,1));
  % Extract profile at a particular y
  [G2,H2,V2,K2,XYP2]=localization2D(XY,[],[],method,[],[],[],...
	    XYP1(XYP1(:,2)==XYP1(wheir,2),:));
  demoplots2(G1,H1,V1,K1,XYP1,XY,G2,H2,V2,K2,XYP2)
elseif strcmp(XY,'demo4')
  [x,y]=blob;
  XY=[x(:) y(:)];
  defval('method','RS')
  % Get the functions on the inscribed data grid
  [G1,H1,V1,K1,XYP1]=localization2D(XY,[],[],method);

  % Now find somewhere specific
  wheir=round(size(XYP1,1)/2);
  wheir=ceil(rand*size(XYP1,1));
  [G2,H2,V2,K2,XYP2]=localization2D(XY,[],[],method,[],[],[],...
	    XYP1(XYP1(:,2)==XYP1(wheir,2),:));
  demoplots2(G1,H1,V1,K1,XYP1,XY,G2,H2,V2,K2,XYP2)
elseif strcmp(XY,'demo5')
 [G1,H1,V1]=localization2D([],[],[],'RS');
 [G2,H2,V2]=localization2D([],[],[],'GL');
 clf
 a=plot(V1,'bo-'); hold on
 b=plot(V2,'rv-');
 legend('RS','GL')
 ylim([-0.1 1.1])
elseif strcmp(XY,'demo6')
 [x,y]=blob;
 XY=[x(:) y(:)];
 [G1,H1,V1]=localization2D(XY,[],[],'RS',[],[],[0.025 0.025]);
 [G2,H2,V2]=localization2D(XY,[],[],'GL',[],[],[0.025 0.025]);
 clf
 a=plot(V1,'bo-'); hold on
 b=plot(V2,'rv-');
 legend('RS','GL')
 ylim([-0.1 1.1])
elseif strcmp(XY,'demo7')
  % For a square patch
  lrtb=[-1 1 1 -1];
  [left,right,top,bottom]=deal(lrtb(1),lrtb(2),lrtb(3),lrtb(4));
  EKS=[left left   right  right left]';
  WAI=[top  bottom bottom top   top ]';
  XY=[EKS(:) WAI(:)];
  % Shannon number
  N=5;
  % Number of tapers required
  J=1;
  % Calculate just one taper of size 129*129
  [G1,H1,V1]=localization2D(XY,N,J,'GL',[],[],[2 2]/129);
  % Then take the interior points of H1 and you have a space-limited
  % taper with the right Shannon number on the right domain
  % Should rewrite to suggest extraction indices
  H1=H1(65:194,65:194);
  ah(1)=subplot(221)
  imagesc(H1); axis image;
  title('H_1')
  ah(2)=subplot(222);
  imagesc(G1); axis image
  title('G_1')
  shrink(ah(1),2,2)
  ah(3)=subplot(223);
  imagesc(decibel(fftshift(abs(fft2(H1,size(G1,1),size(G1,2)))).^2))
  caxis([-50 0]); axis image
  ah(4)=subplot(224);
  imagesc(decibel(fftshift(abs(fft2(G1))).^2)); 
  caxis([-50 0]); axis image
  % Compare with Slepian
  [E,V]=dpss(129,1.15,1); [V^2 V1]
  EE=E(:,1)*E(:,1)';
  imagesc(decibel(fftshift(abs(fft2(EE,size(G1,1),size(G1,2)))).^2))
  caxis([-50 0]); axis image
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demoplots1(G,H,V,K,XYP,XY)
clf
% Watch for centers in plotting!!
for index=1:size(G,3)
  % Get the plot range
  xmins=minmax(XYP(:,1));
  ymins=minmax(XYP(:,2));
  % Note that what comes out of the calculation has min(x) and min(y) in
  % the top left corner of the matrix - not the usual C11, CMN
  % arrangement  - consider using FLIPUD
  % Plot of the bandlimited space-concentrated functions G
  ah(1)=subplot(2,2,1);
  imagefnan([xmins(1) ymins(1)],[xmins(2) ymins(2)],G(:,:,index)); 
  axis image
  hold on; plot(XY(:,1),XY(:,2),'k'); hold off
  xl(1)=ylabel('x');
  yl(1)=ylabel('y');
  % Plot of the spacelimited band-concentrated functions H
  ah(2)=subplot(2,2,2);
    imagefnan([xmins(1) ymins(1)],[xmins(2) ymins(2)],H(:,:,index)); 
  axis image
  hold on; plot(XY(:,1),XY(:,2),'k'); hold off
  xl(2)=ylabel('x');
  yl(2)=ylabel('y');
  % Plot of the power of G
  ah(3)=subplot(2,2,3);
  imagesc(abs(fftshift(fft2(G(:,:,index)))).^2); 
  axis image
  xl(3)=ylabel('k_x');
  yl(3)=ylabel('k_y');
  % Plot of the power of H
  ah(4)=subplot(2,2,4);
  imagesc(abs(fftshift(fft2(H(:,:,index)))).^2); 
  axis image
  xl(4)=ylabel('k_x');
  yl(4)=ylabel('k_y');
  % title(sprintf(...
  %     'j = %i ; V = %8.3f ; min(|G|) = %8.3e ; max(|G|) = %8.3e',...
  %     index,V(index),min(min(abs(G(:,:,index)))),...
  %     max(max((abs(G(:,:,index)))))))
  % Cosmetic adjustments
  fig2print(gcf,'portrait')
  nolabels(ah(2),2); delete(yl(2))
  nolabels(ah(4),2); delete(yl(4))
  longticks(ah)
  pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function demoplots2(G1,H1,V1,K1,XYP1,XY,G2,H2,V2,K2,XYP2)
clf
% Watch for centers in plotting!!

% Show a random eigenfunction from the first four
sind=ceil(rand*4);

% Get the plot range
xmins=minmax(XYP1(:,1));
ymins=minmax(XYP1(:,2));
% Note that what comes out of the calculation has min(x) and min(y) in
% the top left corner of the matrix - not the usual C11, CMN
% arrangement  
% Plot of the bandlimited space-concentrated functions G
ah(1)=subplot(2,2,1);
imagefnan([xmins(1) ymins(1)],[xmins(2) ymins(2)],G1(:,:,sind)); 
axis image
hold on; plot(XY(:,1),XY(:,2),'k'); hold off
xl(1)=ylabel('x');
yl(1)=ylabel('y');
% Plot of the spacelimited band-concentrated functions H
ah(2)=subplot(2,2,2);
imagefnan([xmins(1) ymins(1)],[xmins(2) ymins(2)],H1(:,:,sind)); 
axis image
hold on; plot(XY(:,1),XY(:,2),'k'); hold off
xl(2)=ylabel('x');
yl(2)=ylabel('y');

% Find the profile in the full matrix
profs=XYP1(:,2)==indeks(XYP2(:,2),1);
[i,j]=ind2sub(size(G1(:,:,1)),find(profs));
i=unique(i);

% Plot the profiles - knowing that they are at constant y
ah(3)=subplot(2,2,3);
plot(XYP2(:,1),G2(:,1),'r-'); hold on
plot(XYP1(profs,1),G1(i,j,1),'b+'); 
plot(XYP2(:,1),G2(:,2),'g-'); 
plot(XYP1(profs,1),G1(i,j,2),'bo'); hold off; axis tight
xl(3)=ylabel('x');
yl(3)=ylabel('g(x)');

ah(4)=subplot(2,2,4);
plot(XYP2(:,1),G2(:,3),'r-'); hold on
plot(XYP1(profs,1),G1(i,j,3),'b+'); 
plot(XYP2(:,1),G2(:,4),'g-'); hold on
plot(XYP1(profs,1),G1(i,j,4),'bo'); hold off; axis tight
xl(4)=ylabel('x');
yl(4)=ylabel('g(x)');

% Cosmetic adjustments
fig2print(gcf,'portrait')
nolabels(ah(2),2); delete(yl(2))
nolabels(ah(4),2); delete(yl(4))
longticks(ah)

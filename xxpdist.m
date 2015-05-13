function dxxp=xxpdist(MDX,NDXP,meth)
% dxxp=XXPDIST(MDX,NDXP,meth)
% dxxp=XXPDIST(MDX,[],meth)
%
% Calculates all pairwise distances between the entries of a
% M-dimensional vector of D-dimensional positions and a different
% N-dimensional vector of D-dimensional positions
%
% INPUT:
%
% MDX       An MxD dimensional array with M vectors in D dimensions, 
%           for example, the M positions [x(:) y(:)] in D=2
%           OR: an Mx1 dimensional vector [x(:)] of x-coordinates
%           defining ONE dimension of a grid when also NDXP is thus defined 
% NDXP      An NxD dimensional array with N vectors in D dimensions,
%           for example, the N different positions [xp(:) yp(:)] in D=2
%           ---> If this is empty, the graph is compared to itself
%           OR: an Nx1 dimensional vector [y(:)] of y-coordinates
% meth      1 using something I pieced together [default unless 3 possible]
%           2 my original way of doing things in D=2
%           3 Matlab's built-in way of doing things (needs STATS toolbox)
%           4 as found on http://tinyurl.com/xxpdist
% 
% OUTPUT:
%
% dxxp      An MxN array with the distances between every pair of points
%
% EXAMPLE:
%
% xxpdist('demo1') % Produces no errors if all goes well
%
%% Compute one gridded set of distances from the first corner point
% x=[0:11]*30; y=[0:9]*40;
% [MDX,NDXP]=meshgrid(x(:),y(:));
% dxxp1=xxpdist([MDX(:) NDXP(:)],[MDX(1) NDXP(1)]);
% imagesc(reshape(dxxp1,length(y),length(x))); axis xy
% dxxp7=xxpdist([MDX(:) NDXP(:)],[MDX(7) NDXP(7)]);
%
% SEE ALSO:
% 
% KERNELC2D, LOCALIZATION2D, SSPDIST, SSP21SP
%
% Last modified by fjsimons-at-alum.mit.edu, 01/15/2014

% Not tested for complex input, but the hope is that it works
% For non-integers there might be some rounding error

% Initialize
defval('MDX',[])
defval('NDXP',[])

if ~isstr(MDX)
  if prod(size(MDX))==length(MDX) && prod(size(NDXP))==length(NDXP)
    % This means that both together define one D=2 grid of coordinates
    % Obviously could easily build in D=3 case with ZZ but then again
    [MDX,NDXP]=meshgrid(MDX(:),NDXP(:));
    % Unwrap and pass through again
    MDX=[MDX(:) NDXP(:)];
    NDXP=[];
  end

  % If the STATS toolbox is available, use that one
  if help('stats'); defval('meth',3); end

  % Prepare to do it
  defval('meth',1)

  if ~isempty(NDXP)
    % The numbers of dimensions, D, must be the same of course
    if ~all(size(MDX,2)==size(NDXP,2))
      error('The number of dimensions must be the same in both arrays')
    end
  end

  t=tic;
  % Do it
  switch meth 
   case 1
    % We expand the distance as the square root of the dot product
    % Matlab's DOT function isn't cutting it for us
    % Matlab's BSXFUN function lets us expand the singleton dimensions
    if ~isempty(NDXP)
      dxxp=sqrt(bsxfun(@plus,sum([MDX.*conj(MDX)],2),sum([NDXP.*conj(NDXP)],2)')...
		-2*(MDX*conj(NDXP')));
    else
      dxxp=-2*(MDX*conj(MDX'));
      dxxp=sqrt(dxxp-bsxfun(@plus,diag(dxxp),diag(dxxp)')/2);
    end
   case 2
    if isempty(NDXP)
      % Can't make this any shorter if it's a grid onto itself, so duplicate
      NDXP=MDX;
    end
    if size(MDX,2)~=2 || size(NDXP,2)~=2
      error(sprintf('Method %s only works in 2 dimensions',meth))
    end
    % Obviously could easily build in D=3 case with ZZ but then again
    [XX,XXP]=meshgrid(NDXP(:,1),MDX(:,1));
    [YY,YYP]=meshgrid(NDXP(:,2),MDX(:,2));
    % Could do this with BSXFUN(@minus as in PDIST
    dxxp=sqrt((XX-XXP).^2+(YY-YYP).^2);
   case 3
    if ~isempty(NDXP)
      dxxp=pdist2(MDX,NDXP,'euclidean');
    else
      % The next line takes only about half the time since it does half the work
      dxxp=pdist(MDX,'euclidean');
      % Compare this with TRILOS and TRILOSI... this doubles the time again 
      dxxp=squareform(dxxp);
      % The indexing is, apparently D(I,J)=DD((I-1)*(M-I/2)+J-I)
    end
   case 4
    % Is the DOT function in Matlab somehow more efficient than case 1?? 
    if ~isempty(NDXP)
      dxxp=sqrt(bsxfun(@plus,dot(MDX,MDX,2),dot(NDXP,NDXP,2)')-2*(MDX*conj(NDXP')));
    else
      dxxp=sqrt(bsxfun(@plus,dot(MDX,MDX,2),dot(MDX,MDX,2)')-2*(MDX*conj(MDX')));
    end
   otherwise
    error('Specify a valid method: 1, 2, 3 or 4')
  end

  % Report speed
  disp(sprintf('%s took %f seconds (method %i)',...
               upper(mfilename),toc(t),meth))

  % Provide output
  varns={dxxp};
  varargout=varns(1:nargout);
elseif strcmp(MDX,'demo1')
  x=[0:63]*30; y=[0:127]*40;
  % One set of 2-dimensional coordinates yet to be made
  dxxp11=xxpdist(x,y,1);
  dxxp21=xxpdist(x,y,2);
  dxxp31=xxpdist(x,y,3);
  dxxp41=xxpdist(x,y,4);
  difer(dxxp11-dxxp21,[],[],NaN)
  difer(dxxp11-dxxp31,[],[],NaN)
  difer(dxxp11-dxxp41,[],[],NaN)
  difer(dxxp21-dxxp31,[],[],NaN)
  difer(dxxp31-dxxp41,[],[],NaN)

  [X,Y]=meshgrid(x+1,y+2);
  xp=[0:31]*10; yp=[0:101]*20; [XP,YP]=meshgrid(xp,yp);
  % Two sets of 2-dimensional coordinates fully defined
  dxxp1=xxpdist([X(:) Y(:)],[XP(:) YP(:)],1);
  dxxp2=xxpdist([X(:) Y(:)],[XP(:) YP(:)],2);
  dxxp3=xxpdist([X(:) Y(:)],[XP(:) YP(:)],3);
  dxxp4=xxpdist([X(:) Y(:)],[XP(:) YP(:)],4);
  difer(dxxp1-dxxp2,[],[],NaN)
  difer(dxxp1-dxxp3,[],[],NaN)
  difer(dxxp1-dxxp4,[],[],NaN)
  difer(dxxp2-dxxp3,[],[],NaN)
  difer(dxxp3-dxxp4,[],[],NaN)
  % One set of 2-dimensional coordinates already defined
  dxxp1=xxpdist([X(:) Y(:)],[],1);
  dxxp2=xxpdist([X(:) Y(:)],[],2);
  dxxp3=xxpdist([X(:) Y(:)],[],3);
  dxxp4=xxpdist([X(:) Y(:)],[],4);
  difer(dxxp1-dxxp2,[],[],NaN)
  difer(dxxp1-dxxp3,[],[],NaN)
  difer(dxxp1-dxxp4,[],[],NaN)
  difer(dxxp2-dxxp3,[],[],NaN)
  difer(dxxp3-dxxp4,[],[],NaN)
  % Check that the first set and the last set are identical
  difer(dxxp1-dxxp11,[],[],NaN)
  difer(dxxp2-dxxp21,[],[],NaN)
  difer(dxxp3-dxxp31,[],[],NaN)
  difer(dxxp4-dxxp41,[],[],NaN)
end

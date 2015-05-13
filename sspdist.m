function varargout=sspdist(s,sp,siz,incs)
% [dssp,di,dj,dk]=SSPDIST(s,sp,siz,incs)
% 
% Calculates the distances between matrix entries
%
% INPUT:
%
% s,sp         Running indices into a two-or three-dimensional matrix
%              s   can only be one value, 1x1
%              sp  can be a whole set of values, Mx1 or 1xM
% siz          The SIZE of that matrix (two- or three-d)
% incs         The unit-bearing distance increments (two or three) [defaulted]
%
% OUTPUT:
%
% dssp         The unit-bearing Euclidean distances between s and sp
% di,dj,dk     The row, column, depth steps separating the indices 
%
% EXAMPLE:
%
% sspdist('demo1') % Random entries, no output if working as advertised
% sspdist('demo2',s,fl) % Attempt at reorganizing indices [defaults apply]
%
% SEE ALSO:
%
% XXPDIST, MANDIST, SSP21SP
%
% Last modified by fjsimons-at-alum.mit.edu, 01/15/2014

if ~isstr(s)
  % Default metric is unit, up to three dimensions
  defval('incs',[1 1 1])

  % Correct dimensions, whichever way they came
  sp=sp(:);
  if length(s)~=1
    error('First argument can only be one linear index')
  end
  if length(incs)<2
    error('There must be at least two increment values')
  end

  t=tic;
  % Subscript indices of the first and second linear indices
  % Subscript indices between the first and second linear indices
  if length(siz)==2
    [issp,jssp]=ind2sub(siz,[s ; sp]);
    dk=0;
    incs(3)=0;
  elseif length(siz)==3
    [issp,jssp,kssp]=ind2sub(siz,[s ; sp]);
    dk=kssp(2:end)-kssp(1);
  end
  di=issp(2:end)-issp(1);
  dj=jssp(2:end)-jssp(1);
  % The unitless "Manhattan distances" (in blocks and streets) are now
  mssp=sum(abs([di dj]),2)+abs(dk);
  % The unit-bearing "Manhattan distance" is
  ussp=sum(abs([incs(1)*di incs(2)*dj]),2)+abs(incs(3)*dk);

  % But we really want the unit-bearing Euclidean distances
  % We can't separate out dk, won't use BSXPLUS or REPMAT, so... 
  if dk==0
    dssp=sqrt(sum([(incs(1)*di).^2 (incs(2)*dj).^2],2));
  else
    dssp=sqrt(sum([(incs(1)*di).^2 (incs(2)*dj).^2 (incs(3)*dk).^2],2));
  end

  disp(sprintf('%s took %f seconds',upper(mfilename),toc(t)))

  % Output
  varns={dssp,di,dj,dk};
  varargout=varns(1:nargout);
elseif strcmp(s,'demo1')
  % Make a randomly sized unit grid and check the distance from all grid
  % points to a random entry in it, and compare SSPDIST to XXPDIST
  siz=randij; ri=randi(prod(siz)); [zi,zj]=ind2sub(siz,ri);
  str=sprintf('Distances to element %i (%i,%i) of a %ix%i matrix',...
              ri,zi,zj,siz);
  disp(str)
  a=sspdist(ri,1:prod(siz),siz)';
  b=rindeks(xxpdist(1:siz(2),1:siz(1)),ri);
  % Distances should be the same
  diferm(a,b)
  % Zero should appear
  diferm(a(ri))
  diferm(b(ri))
  clf
  ah(1)=subplot(121);
  imagesc(reshape(a,siz)); % axis image
  ah(2)=subplot(122);
  imagesc(reshape(a,siz)); % axis image
  supertit(ah,str)
  set(ah,'xtick',unique([1 zj siz(2)]),...
         'ytick',unique([1 zi siz(1)]))
  longticks(ah)
elseif strcmp(s,'demo2')
  % Pick out the loop selector (third argument)
  defval('siz',[])
  fl=siz; 
  defval('fl',1)
  
  % Now construct the experiment with a random size
  siz=randij(2^max(3,randi(11))); 
  M=siz(1); N=siz(2);
  % Perhaps change your mind on this
  % M=30; N=20;
  % M=21; N=26;
  % M=480;
  % N=313;
  MN=M*N;
  siz=[M N];
  
  % Extra input - second argument a specific input, third the loop selector
  defval('sp',[])
  s=sp; clear sp
  defval('s',randi(MN))

  % Now tell us what we are about to do
  [zi,zj]=ind2sub(siz,s);
  str=sprintf('Distances to element %i (%i,%i) of a %ix%i matrix',...
              s,zi,zj,siz);
  disp(str)

  % The complete equivalent set of one to all distances
  dsspT=sspdist(s,1:MN,[M N],[1 1]);

  % Determine the indices into the first row of the matrix
  [inds,pm]=ssp21sp(s,[M N],fl);

  if MN<10000
    chex='hard';
  else
    chex='easy';
  end
  
  switch chex
   case 'easy'
    % Easy check (only requires the one to all set)
    dssp1=sspdist(1,1:MN,[M N],[1 1]);
    diferm(dssp1(inds),dsspT)
   case 'hard'
    % Hard check (requires the all-to-all set)
    % The complete equivalent set of all to all distances
    Dij=xxpdist(1:N,1:M);
    % These have to hold 
    diferm(Dij(s,:),dsspT(:)')
    diferm(Dij(1,inds),dsspT(:)')
    diferm(Dij(1,inds),Dij(s,:))
    % Development inspection of the arrays filled so far
    ches=[Dij(1,inds(~~inds))' dsspT(~~inds)]';
    diferm(ches(1,:),ches(2,:))
  end
  
  clf
  spy(reshape(inds,siz))
  shrink(gca,1.1,1.1)
  ts=title(str); movev(ts,-M/10)
  set(gca,'xtick',1:siz(2),'ytick',1:siz(1)); longticks(gca)
  hold on
  plot([zj zj],ylim,'LineS','-','Color',grey)
  plot(xlim,[zi zi],'LineS','-','Color',grey)
  hold off
  drawnow
end

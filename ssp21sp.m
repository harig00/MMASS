function [inds,posmul]=ssp21sp(s,siz,fl)
% [inds,posmul]=SSP21SP(s,siz,fl)
%
% Expresses the distance between a linearly indexed grid node and all the
% other grid nodes as indices into the vector that contains the distances
% between the first linear index of the grid and all the other nodes
%
% INPUT:
%
% s          The linear index into a two-dimensional array
% siz        The size of the array (rows, columns)
% fl         0 Runs a non-redundant full loop with some extra overhead [default]
%            1 Runs a redundant full loop
%            
%            2 Runs a short loop for all positions PAST the given one
%
% OUTPUT:
%
% inds       The indices that map Dij(1,inds)=Dij(s,:) in the distance
%            matrix Dij that contains all pairwise distances between the
%            elements s and sp in the grid
% posmul     The relevant index positions and their multiplicity
%
% SEE ALSO:
%
% SSPDIST, XXPDIST
%
% Last modified by fjsimons-at-alum.mit.edu, 01/15/2014

defval('fl',1)

% Keep time
t=tic;

% The size of the array
M=siz(1);
N=siz(2);
MN=M*N;

% A handful of special cases
if s==1
  inds=1:MN;
  posmul=[inds ; ones(size(inds))];
  return
end

% Allocate index array
inds=zeros(1,MN);
%-[fl==3]*[s+1]);

if s==MN
  inds=MN:-1:1;
  posmul=[inds ; ones(size(inds))];
  return
end

% Only one remaining case after contraction of loops below is when max/min
% don't commute which is when there is a NaN from 0/0, on the last grid row,
% i.e. when a=-1 and d=0 it fails to fill out to the left since we picked
% the max-min order to favor the right side... so prepare for that
if rem(s-1,M)+1==M
  inds(s-M:-M:2)=[1:round([s-M-1]/M)]*M+1;
end

if s>1 && s<MN
  % Row number of s in the grid
  sr=rem(s-1,M)+1;
  % Column number of s in the grid
  sc=ceil(s/M);
  % Symmetry across the main diagonal
  inds(1)=s;
  % The main diagonal
  inds(s)=1;
  % Last column is the first column in reverse
  inds(MN)=MN-s+1;

  % Map all "diagonals" to lower-right quadrant "diagonals"
  redfax=0;
  
  % Get the unique combinations from the below two loop variables
  aa=-sc:N-sc;
  dd=-sr:M-sr;

  % if fl==2 should those be from only 0 in both of them, and then use symmetry

  if fl==0
    % Unique step sizes in their lowest increment and add zero
    [A,D]=meshgrid(aa,dd);
    % The GCD is very expensive
    kpt=gcd(A,D)==1;
    aadd=[A(kpt) D(kpt) ; 0 0];
    % Explicitly non-redundant loop - same result, some overhead
    for ind=1:size(aadd,1)
      a=aadd(ind,1); d=aadd(ind,2);
      js=s+a*M+d:a*M+d:max(2,min(MN-1,[sc+a*(M-sr-[d<0]*[M-1])/d-[d<0]]*M+[d<0]));
      % Any of them we already had?
      redfax=redfax+sum(~~inds(js));
      inds(js)=[1:length(js)]*[abs(a)*M+abs(d)]+1;
    end
  elseif fl==1
    % Explicit full loop - looks the same, but there is redundancy
    for a=aa
      for d=dd
	js=s+a*M+d:a*M+d:max(2,min(MN-1,[sc+a*(M-sr-[d<0]*[M-1])/d-[d<0]]*M+[d<0]));
	% Any of them we already had?
	redfax=redfax+sum(~~inds(js));
	inds(js)=[1:length(js)]*[abs(a)*M+abs(d)]+1;
      end
    end
  end
  % Report on redundancy and check timing
  if redfax>0 && 1==3
    disp(sprintf('%i out of %i are repeated hits, time taken %6.4fs',...
	       redfax,MN,toc(t)))
  end
end

% Check that every element has been hit as it should
diferm(sum(~~inds),MN)

% If second output is requested
if nargout==2
  pos=unique(inds);
  posmul=hist(inds,pos);
  % And the reconstruction must be happening
  diferm(gamini(pos,posmul),sort(inds))
  % Return both the position and its multiplicity
  posmul=[pos; posmul];
end
  

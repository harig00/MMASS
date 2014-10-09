function varargout=surfaz(diro,genx,varargin)
% SURFAZ(diro,genx)
% SURFAZ(diro,genx,nbins)
% SURFAZ(diro,genx,[],1) % Regular grid
% [maxn,minn]=SURFAZ(...)
% demat=SURFAZ(...,1) % Returns full density matrix
%
% Reads the matrices of SURTRACE and calculates azimuthal coverage
%
% INPUT:
%
% diro          A given directory string
% genx          Specifies grid directory inside $GRIDS
% mode          0 one after the other
%               1 all together
% nbins         Number of optional bins to make the color scale
%
% Last modified by fjsimons-at-mit.edu 08/20/2000

[cumnum,gcd,azm,celnr,c11,cmn,lat,dlon,nmr,lonlat,indi,refarea]...
    =lcsg(diro,genx);

if nargin>=4 & varargin{2}==1
  disp('Making regular grid')
  nmr=repmat(max(nmr),length(nmr),1);
  dlon=repmat(min(dlon),length(dlon),1);
end

concel=unique(celnr);
% Such that fillauth(lat,dlon,nmr,c11,concel,1)
% blackens the entire resolved area
binc=[45 135 225 315];

for index=1:length(concel)
  % How many in each quadrant for every concerned cell?
  n(index,:)=hist(azm(celnr==concel(index)),binc);
  % rose(azm(celnr==concel(index))*pi/180,4)
  % pause(0.1)
end

% But the two later bins should be counted with the two former
n=[n(:,1:2)+n(:,3:4)];

% Take minimum of both and use that as criterion
qual=min(n')';
% Not the zeros of course
concel=concel(~~qual);
qual=qual(~~qual);

if nargin>=4 & varargin{2}==1
  demat=full(sparse(1,concel,qual,1,sum(nmr)));
  demat=reshape(demat,nmr(1),length(lat)-1)';
  if nargout
    varargout{1}=demat;
  end
  else
    if nargout
      varargout{1}=min(qual);
      varargout{2}=max(qual);
    end
end

% Remove 5 outliers... this obviously ad hoc...
ab=sort(qual(:));
totgcd(qual>=ab(end-5))=ab(end-6);
disp('5 outliers removed in azimuthal coverage')

if nargin>=3 & ~isempty(varargin{1})
  qual=ceil(scale(qual,[0.01 varargin{1}]));
end

colormap(flipud(gray))
fillauth(lat,dlon,nmr,c11,concel,qual);
plotcont2(c11,cmn)


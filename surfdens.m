function varargout=surfdens(diro,genx,mode,varargin)
% SURFDENS(diro,genx,mode)
% SURFDENS(diro,genx,mode,nbins)
% SURFDENS(diro,genx,mode,nbins)
% SURFDENS(diro,genx,mode,[],1) % Regular grid
% [mindens,maxdens,ch]=SURFDENS(...)
% demat=SURFDENS(...,1) % Returns full density matrix
%
% Reads the matrices of SURTRACE and calculates ray path density in 1/km
%
% INPUT:
%
% diro          A given directory string with the inversion matrices
% genx          Specifies grid directory inside $GRIDS
% mode          0 one after the other
%               1 all together
% nbins         Number of optional pathlength bins
%
% Last modified by fjsimons-at-mit.edu, 08/23/2007

[cumnum,gcd,azm,celnr,c11,cmn,lat,dlon,nmr,lonlat,indi,refarea]...
    =lcsg(diro,genx);

if nargin>=5 & varargin{2}==1
  disp('Making regular grid')
  nmr=repmat(max(nmr),length(nmr),1);
  dlon=repmat(min(dlon),length(dlon),1);
end

switch mode
  case 1
    totgcd=sparse(1,celnr,gcd);

%    Number of rays hitting the cell
%    nrayhit=sparse(1,celnr,ones(size(celnr)));
%    [jk,celnr,nrayhit]=find(nrayhit);
%    colormap(flipud(gray))
%    fillauth(lat,dlon,nmr,c11,celnr,nrayhit);
%    plotcont(c11,cmn)
%
    [jk,celnr,totgcd]=find(totgcd);

    colormap(flipud(gray))
    migcd=min(totgcd);
    magcd=max(totgcd);

    % Remove five outliers
    ab=sort(totgcd(:));
    totgcd(totgcd>=ab(end-5))=ab(end-6);
    disp('5 outliers removed in path density')

    if nargin>=4 & ~isempty(varargin{1})
      totgcd=ceil(scale(totgcd,[0.01 varargin{1}]));
    end    

    fillauth(lat,dlon,nmr,c11,celnr,totgcd);
    [jk,ch]=plotcont(c11,cmn);
  case 0
    counter=1;
    for index=1:length(indi)
      clf
      fillauth(lat,dlon,nmr,c11,...
	  celnr(indi(index,1): indi(index,2)),...
	  azm(indi(index,1): indi(index,2)));
      title(num2str(counter))
      counter=counter+1;
      [jk,ch]=plotcont(c11,cmn);
      pause
    end
end

if nargin>=5 & varargin{2}==1
  demat=full(sparse(jk,celnr,totgcd,1,sum(nmr)));
  demat=reshape(demat,nmr(1),length(lat)-1)';
  if nargout
    varargout{1}=demat;
  end
  else
    if nargout
      varargout{1}=migcd/refarea;
      varargout{2}=magcd/refarea;
      varargout{3}=ch;
    end
end

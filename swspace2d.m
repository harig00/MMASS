function swspace2d(method,scalem)
% SWSPACE2D(method,scalem)
%
% FIGURE 3.1 of SIMONS & WANG, Int. J. Geomath. (2011)
% 
% Spatial rendition of the radial part of the 2D-disk tapers scaled to
% the unit interval. Checks that the finite-length effects are not so
% grave by ensuring that the "small" Fourier components could be "zero"
% after all.
%
% INPUT:
%
% method   'SE' by Slepian "extension" [default]
%          'GL' by direct Gauss-Legendre integration
% scalem   1 Scales solution for weightless orthogonality 
%          0 Scales solution for orthogonality with weight x [default]
%
% SEE ALSO:
%
% SWDISK, SWSPECTRAL2D, SDWSPACE, SDWSPECTRAL
%
% Last modified by fjsimons-at-alum.mit.edu, 09/20/2011

% The method of computation
defval('method','SE');
% The method of scaling, if 0 don't do sqrt(x) scaling of Slepian
defval('scalem',0);

% Number of tapers per order
NM=4;
% Maximum order considered
M=4;
% Shannon number in question
N=42;
% The abscissas - as long as you like to display though see the
% wavenumber-domain representation and SWSPECTRAL2D for the implications
% of changing these values
x=linspace(0,5,2^12);

% This only if you do anything in the wavenumber domain
[fax,selekt]=fftaxis1D(x,length(x),max(x));
kax=2*pi*fax;
% Limiting wavenumber in question
K=2*sqrt(N);
% Bandlimiting array on the full set of wavenumbers
selK=abs([kax ; -flipud(kax(2:end-1))])>K;

% Start the computation and the figure
clf
[ah,ha,H]=krijetem(subnum(M+1,NM));
for m=0:M
  % The above are the unscaled functions "phi" or thus "g"
  [E1,V1{m+1},Nm,c,x,E2]=swdisk(m,N,NM,[],x,method,scalem);
  % This is the integral over the restricted interval
  for inds=1:length(V1{m+1})
    werisit=x(:)<=1 & ~isnan(E1(:,inds));
    difer(2*pi*trapeze(x(werisit),...
		       x(werisit)'.*E1(werisit,inds).^2)-V1{m+1}(inds))
  end
  E1m(m+1)=max(abs(E1(:)));
  if strcmp(method,'GL') && scalem==1
    % Calculate the truly bandlimited spectrum
    FE1=fft(E1); FE1(selK,:)=0;
    % And calculate the truly bandlimited eigenfunctions
    E1F=ifft(FE1);
    % If the input was Hermitian you're all set
    difer(imag(E1F))
    % Turns out in this case you almost get the input back - but only
    % when the normalization is unweighted, if not you'd probably have to
    % put the sqrt(x) in there with the ifft somehow.
  end
  % Plot the space-domain solution
  for ondex=1:NM
    axes(ah(m*NM+ondex))
    pm(m+1,ondex)=plot(x,E1(:,ondex),'-','Color',grey(6));
    hold on
    pg(m+1,ondex)=plot(x(x<=1),E1(x<=1,ondex),'-','Color','k');
    if strcmp(method,'GL') && scalem==1
      % Plot the ifft with the small components removed
      pgg(m+1,ondex)=plot(x,E1F(:,ondex),'-','Color','r');
    end
    drawnow   
  end
end

% Cosmetics
set(ah,'xlim',[0 max(x)],'xtick',[0 1 max(x)],'xgrid','on','ygrid','on',... 
       'ylim',[-1.2 1.2]*max(E1m))
nolabels(ah(1:end-NM),1)
nolabels(ha(M+2:end),2)
longticks(ah,1/2)

% Eigenvalue labels
nf=9;
lox={'ll','ur','lr','ul',...
     'll','ur','lr','ul',...
    'll','ur','lr','ul',...
     'll','ur','lr','ul',...
     'll','ur','lr','ul'};
for m=0:M
  for ondex=1:NM
    axes(ah(m*NM+ondex))
    t{m+1,ondex}=sprintf('%s = %9.6f','\lambda',V1{m+1}(ondex));
    if V1{m+1}(ondex)<0.975
      lox='lr';
    else
      lox='ur';
    end
    [bh(m+1,ondex),th(m+1,ondex)]=...
	boxtex(lox,ah(m*NM+ondex),t{m+1,ondex},nf,[],0.75);
  end  
end
set(th,'FontS',nf-1)

for index=NM*M+1:NM*(M+1)
  axes(ah(index))
  xl(index)=xlabel(sprintf('radius %s=r/R','\xi'));
end

for index=1:M+1
  axes(ha(index))
  yl(index)=ylabel(sprintf('m = %i',index-1));
end

serre(H,1/3,'across')

% So here this goes wrong in the second version
for index=1:NM
  axes(ah(index))
  tlb(index)=title(sprintf('%s = %i','\alpha',index));
end

set([pm(:) ; pg(:)],'LineW',1)

set(ah,'Fonts',nf)

fig2print(gcf,'portrait')

% set(gcf,'Color','w','Inv','off')
shrink(ah,1,1/1.1)

if strcmp(method,'GL') && scalem==1
  set([pgg(:)],'LineW',0.5)
  q=input('Delete red lines? [y] ','s');
  defval('q','y')
  if strcmp(q,'y')
    delete(pgg)
  end
end

figdisp([],[],[],0)

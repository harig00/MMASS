function swspectral2d(method,scalem)
% SWSPECTRAL2D(method,scalem)
%
% FIGURE 3.2 of SIMONS & WANG
%
% Spectral rendition of the 2D-disk tapers scaled on the unit interval.
%
% INPUT:
%
% method   'GL' by direct Gauss-Legendre integration
% scalem   1 Scales the solution for weightless orthogonality [default]
%          0 Scales the solution for orthogonality with weight x
%
% SEE ALSO:
%
% SWSPACE2D, SDWSPACE, SDWSPECTRAL
%
% Last modified by fjsimons-at-alum.mit.edu, 04/14/2009

% Number of tapers per order
NM=4;
% Maximum order considered
M=4;
% Shannon number in question - the higher the better match
N=42;
% Limiting wavenumber in question
K=2*sqrt(N);
% The abscissas - the longer the better the Fourier transform
x=linspace(0,5,2^12);
% The method of computation - the only workable option
method='GL';
% The method of scaling, if 0 don't do sqrt(x) scaling of Slepian
defval('scalem',0);

% The limits on the y-axis
yls=[-85 5]; ytix=[-80 -40 0];
yls=[-72 5]; ytix=[-60 -40 -20 0];

% Put a "physical" scale on the wavenumber space in "regions" or "R"
[fax,selekt]=fftaxis1D(x,length(x),max(x));
kax=2*pi*fax;
% So the wavenumbers should go from the smallest resolvable of
% 2*pi/max(x) R^{-1} to the Nyquist of pi/max(x)*(nfft-1) R^{-1}
% and the "limiting" wavenumber should be corresponding to an area of
% pi in units of R^2 thus with the Shannon number N we get K=2*sqrt(N)
% and the dk should be 2*pi*(nfft-1)/max(x)/nfft

% Start the computation and the figure
clf
[ah,ha,H]=krijetem(subnum(M+1,NM));

for m=0:M
  % The "bandlimited" functions
  [E1,V1{m+1},Nm,c,x,E2]=swdisk(m,N,NM,[],x,method,scalem);
  % The "spacelimited" functions
  E2=E1; E2(x>1,:)=0;
  
  % Perform the Fourier transform to the "spectrum"
  E1=abs(fft(E1)).^2;
  E2=abs(fft(E2)).^2;

  % Plot the spectral-domain solution
  for ondex=1:NM
    axes(ah(m*NM+ondex))
    pm(m+1,ondex)=plot(kax,decibel(E2(selekt,ondex)),'-','Color',grey(6));
    hold on
    % Determine deciBel level based on all of the available power
    dbE1=decibel(E1(selekt,ondex));
    pg(m+1,ondex)=plot(kax(kax<=K),dbE1(kax<=K),...
		       '-','Color','k');
    % As we verify in SWSPACE2D, when we actually set the next component
    % to zero, the spatial-domain eigenfunctions aren't very different at
    % all, which lends credency to the idea that indeed these following
    % small and fast declining values can indeed be thought to be zero
    % after all - coming close to "numerical" bandlimitation
    pgg(m+1,ondex)=plot(kax(kax>K),dbE1(kax>K),...
		       '--','Color','k');
    drawnow
    axis tight
  end
end
% Cosmetics
set(ah,'xlim',[0 max(x)*K],'xtick',[0 K max(x)*K],...
       'xticklabel',[0 K max(x)*K]/K,...
       'xgrid','on','ygrid','on',... 
       'ylim',yls,'ytick',ytix)

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
  xl(index)=xlabel(sprintf('wavenumber %s=k/c','\kappa'));
end

for index=1:M+1
  axes(ha(index))
  yl(index)=ylabel(sprintf('m = %i',index-1));
end

serre(H,1/3,'across')

for index=1:NM
  axes(ah(index))
  tlb(index)=title(sprintf('%s = %i','\alpha',index));
end

set([pm(:) ; pg(:)],'LineW',1)
set([pgg(:)],'LineW',0.5)

set(ah,'Fonts',nf)

fig2print(gcf,'portrait')

% set(gcf,'Color','w','Inv','off')
shrink(ah,1,1/1.1)

if strcmp(method,'GL')
  set([pgg(:)],'LineW',0.5)
  q=input('Delete dashed lines? [y] ','s');
  defval('q','y')
  if strcmp(q,'y')
    delete(pgg)
  end
end

figdisp([],[],[],0)


function Sbar=bluros(S,params,xver)
% Sbar=BLUROS(S,params,xver)
%
% Spectral blurring with periodogram of boxcar. If we're talking about a
% matrix, there are three columns, and there can be an eigenvalue check
%
% INPUT:
%
% S       The spectral matrix with its wavenumbers unwrapped, on a
%         wavenumber grid that is refined from the original
% params  Parameters of this experiment, which is all we need
% xver    1 Do the eigenvalue check [default]
%         0 Don't [but really, should worry about it]
%
% OUTPUT:
%
% Sbar    The blurred spectral matrix, interpolated to original dimension
%
% EXAMPLE:
%
% [Hx,th0,p]=simulosl;
% N=256; h=nan(N-1,1); for index=2:N; p.NyNx=[index index]; 
% Hx=simulosl(th0,p);                       
% h(index-1)=var(Hx)*(2*pi)^2/prod(p.dydx); end
% plot(2:N,h); xlim([1 N]); ylim([0 4e-3]); longticks; shrink(gca);
% t=title(sprintf('SIMULOSL with blurs = 2, var(Hx) versus %s^2 = %g','\sigma',th0(1)))
% movev(t,5/1000); grid on; xlabel('data size Nx = Ny');
% ylabel('var(Hx)*(2pi)^2/(dydx)'); 
% set(gca,'ytick',[0 th0(1) indeks(ylim,2)])
% print('-depsc','/u/fjsimons/EPS/simulosl')
%
% Last modified by fjsimons-at-alum.mit.edu, 10/20/2014

% Set defaults
defval('xver',1)

% New, temporary dimensions
blurs=params.blurs;
NyNx2=blurs*params.NyNx;
% Target dimensions
NyNx=params.NyNx;

% This is the periodogram of the boxcar for the old grid 
% interpolated to the new grid so we don't just hit the nodes
% See http://blinkdagger.com/matlab/matlab-fft-and-zero-padding/
% Should do this once and save it. Now let's not forget that in our
% formalism we force the fft/ifft to be unitary
% If params.blurs were to be 1, we would get a single 1 in the center
Fejk=fftshift(abs(fft2(repmat(1/prod(NyNx)/blurs,NyNx),NyNx2(1),NyNx2(2))).^2);

% Make sure it is unitary and norm-preserving
difer(sum(Fejk(:))-1,[],[],NaN)

% Should use the Claerbout helix! convmtx2 needs more memory
% Actually, should look into FFTW. But also limit to halfplane.
% Which we need to convolve now in two D

%disp(sprintf('BLUROS %i %i %i',blurs,NyNx2(1),NyNx2(2)));

if xver==1
  % The unblurred wavenumbers
  [k,dci,dcn,kx,ky]=knums(params);
  % The blurred wavenumbers
  [k2,kzero,dcn2,kx2,ky2]=knums(params,1);
  % Do the check, don't bother with dcn2 as that's done inside knum2
  difer(k(kzero),[],[],NaN)
else
  % The unblurred wavenumbers
  [~,~,~,kx,ky]=knums(params);
  % The blurred wavenumbers
  [~,~,dcn2,kx2,ky2]=knums(params,1);
end
   
% Check that there is no funny business going on in the interpolation
% This case is going to apply when we start from an odd number and double
if [kx(end)-kx2(end)]>0; kx2(end)=kx(end); end
if [ky(end)-ky2(end)]>0; ky2(end)=ky(end); end

% This case is going to apply when we start from an even number and double
if [kx(1)-kx2(1)]<0; kx2(1)=kx(1); end
if [ky(1)-ky2(1)]<0; ky2(1)=ky(1); end

% Fill to original dimensions
Sbar=nan(prod(NyNx),size(S,2));
% Make sure there are no NaNs in the output
for in=1:size(S,2)
  % Later, consider griddedInterpolant
  Sbar(:,in)=indeks(...
      interp2(kx2(:)',ky2(:),...
	      conv2(Fejk,reshape(S(:,in),NyNx2),'same'),...
	      kx(:)',ky(:)),':');
end

% Check that no extrapolation was demanded, effectively
% but know that griddedInterpolant would have EXTRApolated fine
difer(sum(isnan(Sbar(:))),[],2,NaN)

%% If something is wrong, this gets messed up

% Check the eigenvalues of Sbar for being real and positive
if xver==1 && size(Sbar,2)==3
  ntry=10;
  for kindex=1:ntry
    % Some random wavenumber entry
    kk=max(round(rand*prod(NyNx)),1);
    egos=eig([Sbar(kk,1) Sbar(kk,2) ; Sbar(kk,2) Sbar(kk,3)]);
    if ~isreal(egos) && all(egos>0)
      % disp(sprintf('Some imaginary or nonpositive eigenvalues'))
      cpxity(kindex)=100*mean(abs(imag(egos(:))))./mean(abs(real(egos(:))));
    end
  end
  try
    disp(sprintf(...
	'BLUROS: Maximum im/re percentage out of %i tried is %5.0e%s',...
	ntry,max(cpxity),'%'))
  catch
    %disp(sprintf('BLUROS: No problems found'))
  end
  % Then I would say, if it's too big, don't do it but start over again?
  % Or make the convolution happen on a finer grid?
end



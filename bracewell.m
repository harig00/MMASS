function bracewell
% Illustrates some properties of the Fourier transform
% Calculates the Fourier transform of a signal, and plots its real and
% imaginary parts, and its squared absolute value divided by the sampling
% frequency (the one-sided spectral density).
%
% Last modified by fjsimons-at-alum.mit.edu, 6.10.2005

% Length of the sequence, in physical units, say seconds
tmin=0;
tmax=10*pi;
noil=4;
% Number of time points
nsamp=2^9;
% Number of fft points: equal to nsamp
nfft=nsamp;
% Discrete wavelengths present in the signal
wavl=[0.5 0.2];
% Amplitudes of these components
ampl=[1 1];
% Total time
ttim=tmax-tmin;
% Sampling frequency
Fs=nsamp/ttim;

% Time-domain signal
t=linspace(tmin,tmax,nsamp);
% Construct the time-domain signal
s=zeros(size(t));
% Individual sine harmonic components
for index=1:length(wavl)
  s=s+ampl(index)*sin(t*2*pi/wavl(index));
end
% Add noise
s=s+randn(1,length(t))/noil;
% Frequency-domain signal: the Fourier transform
F=fft(s,nfft);
% The first and the last elements are real, the others complex

% Note that the zero-frequency equals the sum of the data
disp(sprintf('mean(s)=   %8.3f ;      F(1)/N= %8.3f ; diff= %8.3e',... 
	     mean(s),F(1)/length(s),abs(mean(s)-F(1)/length(s))))
% And Parseval's relation that the variance of the data is the variance
disp(sprintf('sqnorm(s)= %8.3f ; sqnorm(F)/N= %8.3f ; diff= %8.3e',...
	     s*s',F*F'/length(s),abs(s*s'-F*F'/length(s))))

% Frequency axis
[f,selekt,l]=fftaxis1D(s,nfft,ttim);
% Select the nonredundant components of the "spectrum"
P=abs(F(selekt)).^2;
% Thus, nfft=nsamp, but get imaginary and real, only half of it is unique,
% thus length(real(F(selekt)))+length(imag(F(selekt)))-2=nsamp again
% But don't forget to double the formerly redundant ones
P=[P(1) 2*P(2:end-1) P(end)];
% 'Cause
if abs(sum(P)-F*F')>1e-8
  error('You are in trouble')
end
% Spectral density
S=P/Fs;

% What is the integral of S over all these frequencies?
ig=trapeze(f,S);
disp(sprintf('sqnorm(s)= %8.3f ;   int(S(f))= %8.3f ; diff= %5.3f',...
	     s*s',ig,abs(s*s'-ig)))
% The accuracy here depends on the integration routine, rather

clf

ah=krijetem(subnum(2,2));

% Plot these various things and annotate them
axes(ah(1))
p=plot(t,real(s));
tl(1)=title('Time Domain, Real');
axis([0 tmax minmax(real(s))*1.1])
grid on
xl(1)=xlabel('t, time (s)');
yl(1)=ylabel('signal amplitude (a)');

axes(ah(2))
plot(f,S)
%tl(2)=title('Spectral density');
axis([minmax(f) minmax(S)*1.1])
%movev(tl(2),range(ylim)/15)
grid on
xl(2)=xlabel('f, frequency (s^{-1})');
yl(2)=ylabel('spectral density (a^2s)');

axes(ah(3))
plot(real(fftshift(F)))
tl(3)=title('Fourier transform');
axis([1 length(F) minmax(real(F))*1.1])
grid on
xl(3)=xlabel('frequency number');
yl(4)=ylabel('real part (a)');
xval=100;
xels=unique([1 xval:xval:nfft nfft]);
xels=[xels(diff(xels)>xval/2) xels(end)];
set(gca,'xtick',xels)

axes(ah(4))
plot(imag(fftshift(F)))
tl(4)=title('Fourier transform');
axis([1 length(F) minmax(imag(F))*1.1])
grid on
yl(4)=ylabel('Imaginary part (a)');
xl(4)=xlabel('frequency number');
set(gca,'xtick',xels)

fig2print(gcf,'portrait')

% Wavelength axis on second plot
xtra=xtraxis1dlin(ah(2),[1/8 1/6 1/4 1/2 1 floor(ttim)],...
		  {num2str(floor(ttim)) '1' '1/2' '1/4' '1/6' '1/8'},1);
axes(xtra)
xlx=xlabel('\lambda, Wavelength (s)');
movev(xlx,-range(ylim)/200)
% Looks bad but prints good
longticks([ah xtra])

figdisp







function [instfreq,t]=installen(x,t,wlen,olap)
% [instfreq,t]=INSTALLEN(x,t,wlen,olap)
%
% Computes the predominant frequency of a signal sampled at regular
% intervals. Following a method suggested by Allen and Kanamori (2003). 
%
% INPUT:
%
% x           The real signal
% t           The times at which the signal is sampled, in s
% wlen        Window length, in seconds (default: 3)
% olap        Window overlap, in percent (default: 95)
%
% OUTPUT:
%
% instfreq    The instantaneous (not angular) frequency
% t           The times at which this is evaluated
%
% EXAMPLE: (Careful with Nyquist frequency and window length!)
%
% t=linspace(0,8*pi,500); f=sin(2*pi*3*t);
% subplot(211); plot(t,f);
% subplot(212); [in,t]=installen(f,t); plot(t,in);
% ylim([0 5])
%
% Last modified by fjsimons-at-alum.mit.edu, 14.04.2006

defval('wlen',3)
defval('olap',95)

% Computes sampling interval
dt=t(2)-t(1);

% Window length, in samples
wlen=round(wlen/dt);

% Overlap, in samples
olap=floor(olap/100*wlen);

% It's a boxcar window in this case
dwin=ones(wlen,1);

% Determine window parameters and truncation of data set
% Number of points in the signal
npts=length(x);
% Effective length of the window
eflen=wlen-olap;
% Check this works
checkit=(npts-olap)/eflen;
nwin=floor(checkit);
disp(sprintf('Number of overlapping data segments: %i',nwin))   
if nwin~=checkit
  disp(sprintf(...
      'Number of  segments is not  integer: %i / %i points truncated',...
      npts-(nwin*eflen+olap),npts))
end
if nwin<=0; error('Data sequence not long enough'); end
% First, divide the signal in overlapping windowed segments
rows=[1:wlen];
cols=[0:(nwin-1)]*eflen;
xsd=detrend(...
    x(repmat(rows(:),1,nwin)+...
      repmat(cols,wlen,1)));
% This is the overlapping windowed data set
% Do not normalize this window... this is not the spectrogram
xsdw=xsd.*repmat(dwin,1,nwin);

% Now do algorithm from Allen and Kanamori
% Don't do the cumlative sum, just return the result at the end of the  
% window. Note that this is unstable due to the division by a small number!
in=sum([repmat(0,1,nwin) ; diff(xsdw,1)/dt].^2,1);
xsdw=sum(xsdw.^2,1);
in=sqrt(in(:)./xsdw(:));

instfreq=in/2/pi;
t=cols(:)*dt;




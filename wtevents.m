function newo=wtevents(findex,ddir,tipe,nvm,nd,pph,intel,thresh)
% WTEVENTS(findex,ddir,tipe,nvm,nd,pph,intel,thresh)
%
% Makes diagnostic plots of a single "seismogram"
% using WAVELET ANALYSIS:
% Add option for filtering.
% Add thresholding in the title.
%
% Plot UL: Time-domain representation of the signal
% Plot UR: Successive approximation of the signal
% Plot LL: Scalogram: power of wavelet coefficients
% Plot LR: Successive reproducibility criteria
%
% INPUT:
%
% findex           [scalar] a running number identifying the file in ddir
%                  [string] a filename string
% ddir             is the directory (with the events identified by STALTA)
%                  [default:~/MERMAID/SACDATA/EVENTS]
% tipe             wavelet type ['CDF']
% nvm              number of vanishing moments [[2 4]]
% nd               number of cascades
% pph              1 Time-domain full rate
%                  2 Time-domain polyphase
%                  3 Z-domain polyphase
%                  4 Lifting implementation
%
% EXAMPLE:
%
% wtevents('Okal.sac',pwd,[],[],[],[],[],1)
% wtevents('Okal.sac',pwd,[],[],[],4,[],0)
% wtevents('BBR.BHZ.sac.t0','/home/fjsimons/EALARMS/20010213.030435/BHZdata')
%
% See also EVENTS, SIGNALS, WTEVENTS2 (and higher numbers)
%
% Last modified by fjsimons-at-alum.mit.edu, 16.11.2005

% Establish defaults
defval('findex',3)
defval('tipe','CDF')
defval('nvm',[2 4])
defval('nd',5)
defval('pph',4)
defval('intel',0)
defval('thresh',1)
defval('ddir','/home/fjsimons/MERMAID/08092007/EVENTS/')

tph{1}='Time-Domain Full Rate';
tph{2}='Time-Domain Polyphase';
tph{3}='Z-Domain Polyphase';
tph{4}='Lifting Implementation';

tph=tph{pph};

% Get data
if ~isstr(findex)
  files=ls2cell(ddir);
  filen=files{findex};
else
  filen=findex;
end

titl=nounder(sprintf('%s / %s %i,%i / %s',filen,tipe,nvm,tph));

clf
ah=krijetem(subnum(2,2));

% TIME-DOMAIN REPRESENTATION OF SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(1))
filenam1=fullfile(ddir,filen);
[x1,h1,t,p]=readsac(filenam1,0,'l');
cof=2;
%sig=x1;
sig=lowpass(x1,1/h1.DELTA,cof,2,2,'butter');
x1=sig;
timax=linspace(h1.B,h1.E,h1.NPTS);
p(1)=plot(timax,sig);
%xll(1)=xlabel(sprintf('%s ; Low-pass 2 Hz ;  %i s selected',...
%		      'Time (s)',ceil(h1.T1-h1.T0)));
xll(1)=xlabel(sprintf('%s ; %i s selected',...
		      'Time (s)',ceil(h1.T1-h1.T0)));
axis tight
yl(1)=ylabel(sprintf('Filtered Amplitude, %3.1f Hz',cof));
hold on
grid on
yli=ylim;
pt01(:,1)=plot(repmat([h1.T0 h1.T1],2,1),[yli(1) yli(2)],'k--');

% SUCCESSIVE APPROXIMATION OF SIGNAL BY WT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ochk=1*mod(length(x1),2);
[a,d,an,dn,xs]=succapp(tipe,nvm,nd,x1(1:end-ochk),pph,intel);
axes(ah(2))
timax=linspace(h1.B,h1.E-h1.DELTA*ochk,h1.NPTS-ochk);
if thresh==1
  [dt,dnz]=threshold(d,dn,[],2);
  [xr,ts]=iwt(a,dt,an,dn,tipe,nvm,pph,intel);
  xrs=[xr{:}];
  [psa,yl]=succaplot(timax,xrs);
else
  [psa,yl]=succaplot(timax,xs);
end
xlim([h1.B h1.E])
yll(2)=ylabel('Minimum Coarseness Included');
xll(2)=xlabel('Time (s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SCALOGRAM
axes(ah(3))
defval('meths',1);

if thresh==1
  [pdy,stp,rgp]=dyadplot(x1(1:end-ochk),a,dt,an,dn,meths,[h1.B h1.E]);
else
  [pdy,stp,rgp]=dyadplot(x1(1:end-ochk),a,d,an,dn,meths,[h1.B h1.E]);
end
hold on
yli=ylim;
pt01(:,2)=plot(repmat([h1.T0 h1.T1],2,1),[yli(1) yli(2)],'k--');
if meths==1
  colormap(flipud(gray(128)))
  caxis([0 3*stp])
end

yll(3)=ylabel('Resolution level (Coarseness)');
xll(3)=xlabel('Time (s)');
set(ah(2:3),'Xgrid','on','Ygrid','on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT COMPRESSION EFFICIENCY
axes(ah(4))
if thresh==1
  [cf,rmse,acor,ener]=compres(x1(1:end-ochk),xrs,an,dnz);
else
  [cf,rmse,acor,ener]=compres(x1(1:end-ochk),xs,an,dn);
end
pr=plot(log(cf),rmse,'-o'); hold on
pa=plot(log(cf),acor,'-s');
pe=plot(log(cf),ener,'-v');
hold off
ylim([0 101])
xlim([floor(indeks(log(cf),1)) ceil(indeks(log(cf),'end'))])
xll(4)=xlabel('ln(# coefficients / # data)');
yll(4)=ylabel('Reproducibility criteria (%)');
ax=xtraxis(ah(4),log(cf),floor(cf),'# coefficients / # data (%)');
set(ah(4),'Xgrid','on','Ygrid','on')
st=supertit(ah(1:2),titl);
movev(st,.3)
longticks([ah ax])
fig2print(gcf,'landscape')
set(pr,'MarkerF','b','MarkerE','b','Color','b')
set(pa,'MarkerF','w','MarkerE','k','Color','k')
set(pe,'MarkerF','r','MarkerE','r','Color','r')
axes(ah(4))
wl=4; if rmse(end-1)<40; wl=2; end
l=legend('RMSE','XCORR','ENER',wl);
axes(ax)

id

if pph==3
  if ~isstr(findex)
    figdisp(sprintf('wtevents_%3.3i',findex))
  else
    figdisp(sprintf('wtevents_%s_l',pref(findex)))
  end
elseif pph==4
  if ~isstr(findex)
    figdisp(sprintf('wtevents_%3.3i_l',findex))
  else
    figdisp(sprintf('wtevents_%s_l',pref(findex)))
  end
end

disp('PRINT WITH ZBUFFER ON')


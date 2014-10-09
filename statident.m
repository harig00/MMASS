function [statio,varargout]=statident(lonlat,tol)
% statio=STATIDENT(lonlat,tol)
% [statio,lon,lat]=STATIDENT(lonlat,tol)
%
% Identifies stations from a database by their name or by their
% coordinates; returns one stationname per location OR returns one
% location per stationname. Returns either exact matches or those
% stations that are within some tolerance of the single target pair.
%
% The data base, stations.mat, as read by the program STATION,  is a
% single structure array, formatted as: 
%
%    sc09: [136.2953 -22.3969]
%    sc08: [131.3672 -22.2748]
%
% INPUT:
%
% lonlat     [lon lat]  Longitude(s) and latitude(s) (decimal degrees) OR
%            {'lonlat'}  Station name(s) (cell)
% tol        Tolerance (km) on the distance to the coordinates 
%
% OUTPUT:
%
% statio     The name of the station OR
%            a matrix with lon and lat
%
% See also STATION, STRUCT2ASCI, LOCIDENT
%
% EXAMPLE: 
%
% statident([241.714  34.019; 241.829  34.148 ])
% statident({'USC' 'PASA'})
% statident(statident({'USC' 'PASA'}))
% [statio,lon,lat]=statident([241.7 33.9],20);
% statident([136 -22],200)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/02/2012

% Get the station list
[s,sn,sc]=station;

% Default is no tolerance
defval('tol',[])

if ~iscell(lonlat) & ~isstr(lonlat)
  if size(lonlat,1)==1 & ~isempty(tol)
    gcdkm=grcdist(lonlat,sc);
    pos=1:length(sn);
    statio=sn(pos(gcdkm<tol));
    varargout{1}=sc(pos(gcdkm<tol),1);
    varargout{2}=sc(pos(gcdkm<tol),2);
  else
    stesj=[];
    for index=1:length(sn)
      stesj=str2mat(stesj,sn{index});
    end
    stesj=stesj(2:end,:);
    statio=[];
    for index=1:size(lonlat,1)
      indes=find(sc(:,1)==lonlat(index,1)  & sc(:,2)==lonlat(index,2));
      if ~isempty(indes)
	statio=str2mat(statio,stesj(indes(1),:));
      else
	statio=str2mat(statio,[]);
      end
    end
    if ~isempty(indes)
      statio=statio(2:end,:);
    end
    if size(statio,1) ~= size(lonlat,1)
      warning('Not all stations found')
    end
  end
else
  if isstr(lonlat)
    %statio=eval(sprintf('s.%s',lonlat));
    statio=s.(lonlat));
  else
    for index=1:length(lonlat)
      %statio(index,:)=eval(sprintf('s.%s',lonlat{index}));
      statio(index,:)=s.(lonlat{index});
    end
  end
end

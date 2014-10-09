function [s,sn,sc]=station
% [s,sn,sc]=STATION
%
% Comes up with a list of stations, names and coordinates
% from a binary Matlab data base file
%
% INPUT:
%
% Nothing.
%
% OUTPUT:
%
% s      A structure array with station names and coordinated
% sn     A cell array with station names
% sc     A matrix with longitude, latitude in decimal degrees
%
% See also STATIDENT, STRUCT2ASCI
%
% Last modified by fjsimons-at-alum.mit.edu, 07/15/2010

% You need an environmental variable, the directory, and the named file
load(fullfile(getenv('IFILES'),'STATIONS','stations'))

% Now do it
stt=struct2cell(stations);
stt=[stt{:}];
stt=reshape(stt,2,length(stt)/2)';

sn=fieldnames(stations);
s=stations;
sc=stt;


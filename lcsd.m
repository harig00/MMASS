function [L,C,S,D]=lcsd(diro)
% [L,C,S,D]=LCSD(diro)
%
% Loads L, C, S and D
L=loadb(fullfile(diro,'L.bin'),'int32');
C=loadb(fullfile(diro,'C.bin'),'int32');
S=loadb(fullfile(diro,'S.bin'),'float32');
D=loadb(fullfile(diro,'D.bin'),'float32');


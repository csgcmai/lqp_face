function [feat,nrows,colsize]=readFaceFeatureFile(fname)
% function reads the binary feature files written by cpp code.
% input --> fname: binary file name
% output <-- feat: features as a row vector
%	     nrows: number of horizontal cells.
%	     ncols: number of vertical cell x codebook size
fid=fopen(fname,'rb');
version=fread(fid,1,'integer*4');
nrows=fread(fid,1,'integer*4');
colsize=fread(fid,1,'integer*4');
feat=fread(fid,nrows*colsize,'real*4');
fclose(fid);

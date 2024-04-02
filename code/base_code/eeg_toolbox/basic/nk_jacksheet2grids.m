function [grids grid_names] = nk_jacksheet2grids(jackFile)
% DESCRIPTION:
%  outputs the grid cell array that is fed into reref.m for data 
%  split with nk_split.m  
%  
% function [grids grid_names] = nk_jacksheet2grids(jackFile)
%
% inputs:
%  (1) jackFile: '/data/eeg/TJ017/docs/jacksheet.txt'  
%
% output:
%  (1) grids: the cell array to feed to reref.m 
%  (2) grid_mames: the names that each cell from 'grids"
%                  corresponds to
%
  
fid = fopen(jackFile);
foo = textscan(fid,'%s%s');
fclose(fid)

channels = foo{1};
names    = foo{2};
stripNames = {};

% get the unique stripNames
for k=1:length(foo{2})
  thisName = names{k};
  stripNames{k} = thisName(regexp(thisName,'[^0-9]'));
end
unStripNames = unique(stripNames);

% make the grids
grids = {};
grid_names = {};
count = 0;
for k=1:length(unStripNames)
  thisStripName = unStripNames{k};
  if strcmp(upper(thisStripName),'EKG');
    continue
  end
  count = count + 1;
  grids{count,1}=find(strcmp(thisStripName,stripNames));
  grid_names{count,1}=thisStripName;
end

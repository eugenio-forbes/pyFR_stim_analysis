function [v,f,eNames] = read_surf_wrapper(varargin)
% [v,f,eNames] = read_surf_wrapper(path2surf)
%   This function returns surface data in a manner that can be plotted by
%   plotsurf and plotsurf_w_elecs_wrapper. Automattically 1 indexes it so
%   no need to do that next. Also returns eNames if you are loading a .mat
%   surf with electrode information (see writesurf_w_elecs)
%
% INPUTS:
%   path2surf -         The path to the surface file for this subject
%                       (lh.pial)
%   ---OR---
%   subj      -         The subject whose brain is to be plotted
%   surfType  -         The type of the surface to be plotted 
%                       (lh.pial, rh.pial, lh.pial-outer-smoothed, etc)
% OUTPUTS:
%   v - vertices (xyz coords)
%   f - faces (representation of patches used for the surface);
%       if there is a fourth column, it indicates the electrode that it
%       has been assigned to
% eNames - electrode names associated with electrodes on the surfae 
%          (e.g., '1-2'; see talStruct.eNames)
% created 11/13 by Ashwin Ramayya (ashwinramayya@gmail.com)
% edited to allow second set of arguments 6/14 by Isaac Pedisich
%                                  (iped@sas.upenn.edu)
if nargin<1 || nargin>2
    error('read_surf_wrapper:improperArguments','read_surf_wrapper requires 1 or 2 arguments');
end

if nargin==1
    path2surf=varargin{1};
    if strcmp(path2surf(end-3:end),'.mat')
        %load surface file with electrode info saved by writesurf_w_elecs
        load(path2surf)
        
        % already 1 indexed
    else
        %load freesurfer surface file
        % load surf
        [v,f] = read_surf(path2surf);
        
        % make it 1 indexed
        f = f+1;
        
        eNames = {};
        
    end
else
    [v,f] = read_surf(fullfile('/data/eeg/freesurfer/subjects',varargin{1},'surf',varargin{2}));
    f = f+1;
end
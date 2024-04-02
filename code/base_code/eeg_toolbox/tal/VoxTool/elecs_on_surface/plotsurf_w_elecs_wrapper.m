function [p] = plotsurf_w_elecs_wrapper(path2surf_w_elecs, eNames, kolList,cView,projOnSurfFlag,elecLocs,mSize)
% [p] = plotsurf_w_elecs_wrapper(path2surf_w_elecs, eNames, kolList,cView,projOnSurfFlag,elecLocs,mSize)
%This function can be used to plot electrodes on the cortical surface of the brain
%plotsurf is a function in the iso2mesh package and has been edited to plot
%subregions.
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc
%
%INPUTS
% path2surf_w_elecs             path to .mat surf file which has been saved
%                               by writesurf_w_elecs
% eNames                        names of electrodes that should be colored
%                               (eg, '1-2', see talStruct.eNames)
%  kolList                      [r g b] associated with each electrode in eNames 
%                               -normalized; e.g. 1 0 0 = red
% cView (optional)             - camra view [az el] (default = [90 0])

% projOnSurfFlag (optional)      if set to 1, it will project electrodes onto the surface                            
                                %if 0, it will display the electrodes
                                %'overlayed' on the surface
% elecLocs (optional)           % matrix of electrode locations
                                %-should be in surfaceRAS (snapped or not)
                                %-only needed if ~projOnSurfFlag
% mSize (optional)           % size of electrodes
                                % -only if ~projOnSurf

% returns
%   p : the handle to the plot created

%parse inputs
if ~exist('projOnSurfFlag','var') || isempty(projOnSurfFlag)
    projOnSurfFlag = 1;
end
if ~exist('cView','var') || isempty(cView)
    cView = [-90 0];
end
if ~exist('mSize','var') || isempty(mSize)
    mSize = 10;
end

% read surface
[v,f,eNames_on_surf] = read_surf_wrapper(path2surf_w_elecs);

% initialize kolList_on_surf
kolList_on_surf = nan(length(eNames_on_surf),3);

% update kolList_on_surf based on queried electrodes and requested colors
% ( only if projElecsOnSurf == 1)
if projOnSurfFlag
    for e = 1:length(eNames)
        surfIdx = strcmp(eNames{e},eNames_on_surf);
        kolList_on_surf(surfIdx,:) = kolList(e,:);
    end
end
% pad  kolList_w_surf with nanes as the first rows;
% this is because f(:,4) == 0 reflects cortical surface 
% that is not associated with any electrode and should not be colored (see
% line 124
kolList_on_surf = [nan nan nan; kolList_on_surf];

p = plotsurf(v, f, kolList_on_surf,'facecolor','white','linestyle','none');

% run through and plot electrodes on the surface
if ~projOnSurfFlag 
    hold all
    for e = 1:size(elecLocs,1)
       plot3(elecLocs(e,1),elecLocs(e,2),elecLocs(e,3),'o','markeredgecolor','k','markerfacecolor',kolList(e,:),'markersize',mSize);
    end
end
view(cView)
camlight;
    
    
    

%end

function hm=plotsurf(node,face,kolList,varargin)
% NOTE: This is edited from the original version to allow for electrode
% plotting -AGR
% hm=plotsurf(node,face,opt)
%
% plot 3D surface meshes
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: node coordinates, dimension (nn,3); if node has a 
%            4th column, it will be used to set the color at each node.
%      face: triangular surface face list; if face has a 4th column,
%            it will be used to separate the surface into 
%            sub-surfaces and display them in different colors;
%            face can be a cell array, each element of the array represents
%            a polyhedral facet of the mesh, if an element is an array with
%            two array subelements, the first one is the node index, the
%            second one is a scalar as the group id of the facet.
%      opt:  additional options for the plotting, see plotmesh
%
% output:
%   hm: handle or handles (vector) to the plotted surfaces
%
% example:
%
%   h=plotsurf(node,face);
%   h=plotsurf(node,face,'facecolor','r');
% 
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
rngstate = rand ('state');
if(nargin>=2)
  randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"
  if(isoctavemesh) randseed=randseed+3; end
  if(~isempty(getvarfrom({'caller','base'},'ISO2MESH_RANDSEED')))
        randseed=getvarfrom({'caller','base'},'ISO2MESH_RANDSEED');
  end
  rand('state',randseed);

  if(iscell(face))
    sc=sparse(10,3); % face colormap
    sc(1:10,:)=rand(3,10)';
    len=length(face);
    newsurf=cell(1);
    % reorganizing each labeled surface into a new cell
    for i=1:len
        fc=face{i};
        if(iscell(fc) && length(fc)>=2)
            if(fc{2}+1>10)
                sc(fc{2}+1,:)=rand(1,3);
            end
            if(fc{2}+1>length(newsurf))
                newsurf{fc{2}+1}={};
            end
            newsurf{fc{2}+1}{end+1}=fc{1};
        else % unlabeled facet is tagged by 0
            if(iscell(fc))
                newsurf{1}{end+1}=cell2mat(fc);
            else
                newsurf{1}{end+1}=fc;
            end
        end
    end
    hold on;
h=[];
    newlen=length(newsurf);

    for i=1:newlen
        if(isempty(newsurf{i})); continue; end
        try 
            subface=cell2mat(newsurf{i}')';
            if(size(subface,1)>1 && ndims(subface)==2)
               subface=subface';
            end
            h=[h patch('Vertices',node,'Faces',subface,'facecolor',sc(i,:),varargin{:})];
        catch
            for j=1:length(newsurf{i})
                h=[h patch('Vertices',node,'Faces',newsurf{i}{j},'facecolor',sc(i,:),varargin{:})];
            end
        end
    end
  else
    if(size(face,2)==4)
        tag=face(:,4);
        types=unique(tag);
        hold on;
        h=[];
    for i=1:length(types)
            if(size(node,2)==3)
                % AGR: edited to allow for kolList inputs
                if sum(isnan(kolList(i,:))) == 3
                     h=[h plotasurf(node,face(find(tag==types(i)),1:3),varargin{:})];
                else
                    h=[h plotasurf(node,face(find(tag==types(i)),1:3),'facecolor',kolList(i,:),'linestyle','none')];
                end
                %h=[h plotasurf(node,face(find(tag==types(i)),1:3),'facecolor',rand(3,1),varargin{:})];
            else
                h=[h plotasurf(node,face(find(tag==types(i)),1:3),varargin{:})];
            end
        end
    else
        h=plotasurf(node,face,varargin{:});
    end
  end
end    
if(~isempty(h)) 
  axis equal;
  if(all(get(gca,'view')==[0 90]))
      view(3);
  end
end
if(~isempty(h) && nargout>=1)
  hm=h;
end

rand ('state',rngstate);

%-------------------------------------------------------------------------
function hh=plotasurf(node,face,varargin)
isoct=isoctavemesh;
if(size(node,2)==4)
if(isoct && ~exist('trisurf','file'))
    h=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),node(:,4),'edgecolor','k',varargin{:});
else
    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),node(:,4),varargin{:});
end
else
if(isoct && ~exist('trisurf','file'))
    h=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'edgecolor','k',varargin{:});
else
    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),varargin{:});
end
end
if(exist('h','var')) hh=h; end

function [isoctave verinfo]=isoctavemesh
%
% [isoctave verinfo]=isoctavemesh
%
% determine whether the code is running in octave
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% output:
%   isoctave: 1 if in octave, otherwise 0
%   verinfo: a string, showing the version of octave (OCTAVE_VERSION)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
verinfo='';
%isoctave=(exist('OCTAVE_VERSION')~=0);
isoctave = 0;
if(nargout==2 && isoctave)
    verinfo=OCTAVE_VERSION;
end
function p=getvarfrom(ws,name)
%
% p=getvarfrom(ws,name)
%
% get variable value by name from specified work-space
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%    ws: name of the work-space, for example, 'base'
%    name: name string of the variable
%
% output:
%    p: the value of the specified variable, if the variable does not
%       exist, return empty array
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

wsname=ws;
if(~iscell(ws))
   wsname=cell(1);
   wsname{1}=ws;
end

p=[];
for i=1:length(wsname)
    isdefined=evalin(wsname{i},['exist(''' name ''')']);
    if(isdefined==1)
        p=evalin(wsname{i},name);
        break;
    end
end

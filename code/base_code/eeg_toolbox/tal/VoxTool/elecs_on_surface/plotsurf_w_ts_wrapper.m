function [p] = plotsurf_w_ts_wrapper(v,f,face_ts,cRange,noElecsColor,nonSigColor)
% function [p] = plotsurf_w_ts_wrapper(v,f,face_ts,cRange)
%This function plots the surface and colors each face based on t-values
%provided in face_ts. cRange sets the saturation. Can plot mean power or
%counts with this function, but must adjust cRange accordingly

%plotsurf is a function in the iso2mesh package and has been edited to plot
%subregions.
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc
%
%INPUTS
% v .... xyz coords of each vertex (from readsurf_wrapper, NOT readsurf)
% f .... indices of vertices from each face (from readsurf_wrapper, NOT
         % readsurf)
%face_ts .... vector of t-stats associated with each face (see
         %ec_plotTsOnAvgSurf)
%cRange  .... color saturation for the plot
%noElecsColor (optional) .... how to color regions without coverage
                         %(default = white, which shows the cortical surface)
%nonSigColor (optional) .... how to color regions that are not-significance
                         %(default = .5 .5 .5).
                         


% returns
%   p : the handle to the plot created

if ~exist('noElecsColor','var')||isempty(noElecsColor)
   noElecsColor = [1 1 1]; 
end
if ~exist('nonSigColor','var')||isempty(nonSigColor)
   nonSigColor = [.5 .5 .5]; 
end

% Identify unique ts and assign a color to each t-stat
u_ts = unique(face_ts(~isnan(face_ts)));
if isempty(u_ts)
    u_kolList = [];
else
    u_kolList = stat2kol(u_ts,cRange,100);
    % t-stats with 0 should show darker cortical surface ([1 1 1])
    if sum(u_ts==0)>0
        u_kolList(u_ts==0,:) = nonSigColor; 
    end
end


% add fourth column to faces so as to tag similarly colored faces
f_tag = zeros(size(face_ts));

for k = 1:length(u_ts)
   idx = face_ts == u_ts(k);
   f_tag(idx) = k;  
end

f = [f, f_tag'];


% if (f:,4) == 0, it should color the surface white/black
%as it is not associated with a t-stat (see
% line 124
u_kolList = [nan nan nan; u_kolList];

p = plotsurf(v, f, u_kolList,'facecolor',noElecsColor,'linestyle','none');

%view;camlight    
  
%remove glossiness
 h=findobj('type','patch');
 for k=1:length(h);
    set(h(k),'SpecularStrength',.1, ...
	  'DiffuseStrength',.6, ...
	  'SpecularColorReflectance',0, ...
	  'AmbientStrength',.45);
end   

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

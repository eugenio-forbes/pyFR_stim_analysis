function colorBrainSurface(elecLocs, surfL, surfR, radius, colors)
%COLORBRAINSURFACE 
%This function plots a surface, and colors regions nearby electrode locations
%Inputs
%elecLocs           %matrix of xyz coordinates associated with electrodes
                    %(in surfaceRAS space) 
%surfL (optional)   %path to left surface (pial)
%surfR (optional)   %path to right surface (pial)
%radius             %radius to use when identifying regions of the surface
                    %near the electrode
%colors             %colors associated with each electrode 

%Written by Isaac 11/22/2013. Preamble by ashwin ramayya
                    
hold all
if exist('surfR','var') && ~isempty(surfR)
    [v_R,f_R] = read_surf(surfR);
    plotsurf_wrapper(v_R, f_R+1);
end
if exist('surfL','var') && ~isempty(surfL)
    [v_L,f_L] = read_surf(surfL);
    plotsurf_wrapper(v_L, f_L+1);
end

for i=1:size(elecLocs,1)
    fprintf('%d ',i)
    point = elecLocs(i,:);
    if exist('surfR','var') && ~isempty(surfR)
        indices_R = 1:length(v_R);
        pointRep = repmat(point, size(v_R,1),1);
        distances = sqrt(sum((pointRep-v_R).^2,2));
        relevant_indices = indices_R(distances<radius);
        f_R_relevant = f_R(any(ismember(f_R, relevant_indices),2),:);
        if size(f_R_relevant,1)>0
            hold all
            plotsurf_wrapper(v_R, f_R_relevant+1, colors(i,:));
        end
    end
    
    if exist('surfL','var') && ~isempty(surfL)
        pointRep = repmat(point, size(v_L,1),1);
        distances = sqrt(sum((pointRep-v_L).^2,2));
        indices_L = 1:length(v_L);
        relevant_indices = indices_L(distances<radius);
        f_L_relevant = f_L(any(ismember(f_L, relevant_indices),2),:);
        if size(f_L_relevant,1)>0
            hold all
            plotsurf_wrapper(v_L, f_L_relevant+1, colors(i,:));
        end
    end
    drawnow;
    
end



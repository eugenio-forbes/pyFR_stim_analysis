function [elecNums, noLocFound, distToHippHead, js_error, hipp_error]=get_mtl_leads(subj,region,neuroRadFlag)
% DESCRIPTION:
%  This function returns the electrode numbers associated with 
%
% FUNCTION
%    [elecNums, noLocFound]=get_mtl_leads(subj,region,neuroRadFlag)
%
% INPUT:
%   subj ..... the subject you want (e.g., 'TJ069')
%   region (optional) ... mtl region of interest
%            'hipp' = hippocampal leads (default)
%             
%   neuroRadFlag (optional) .... which localization scheme to use
%                           0 (default)..... joel stein preferred
%                            tries neurorad_localization.txt, but get
%                            hippo_only.txt otherwise
%                           1 ..... joel stein only  
%                           2 ..... hippleads only
% OUTPUT:
%  elecNums ... from the jacksheet (correspond to 'eNames' in subjTalStruct)
%  noLocFound ... is true if unable to find the localization infor
                %  requested. Use taliarach daemon infor otherwise
% distToHippHead...distance between hippocampus head and electrode along long axis of hippocampus
% js_error  .... joel stein loc not found
% hipp_error ... hipp_leads.info not found

% NOTE:  
%
% Written by ashwin g. ramayya (ashwinramayya@gmail.com)
if ~exist('region','var') || isempty (region)
    region = 'hipp';
end
if ~exist('neuroRadFlag','var') || isempty (neuroRadFlag)
    neuroRadFlag = 0;
end
switch upper(region)
    case{'HIPP','H'} 
        [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,'HIPP',[]);
        [hipp_elecNums,hipp_error] = gethippleads(subj,0);
    
    otherwise
        %EC, PHC, PRC, FUSIFORM, ITG, MTG, STG 
        [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,upper(region),[]);
        hipp_elecNums =  []; hipp_error = 1;
end
   
% parse neuroRad flag
if neuroRadFlag == 2 
    elecNums = hipp_elecNums;
    noLocFound = hipp_error;
else
    % get js_enums
    elecNums = js_elecNums;
    noLocFound = js_error;
    
    if js_error && (neuroRadFlag == 0)
       elecNums = hipp_elecNums;
       noLocFound = hipp_error;
    end
end


% case{'PHC','parahippocampalcortex'} 
% case{'EC','entorhinalcortex'} 
% [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,'EC',[]);
%  hipp_elecNums =  []; hipp_error = 1;
% case{'PRC','perirhinalcortex'} 
% [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,'PRC',[]);
% hipp_elecNums =  []; hipp_error = 1;
% case{'MTL'}
% [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,'MTL',[]);
% case{'FUSIFORM'}
% [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,'FUSIFORM',[]);
% case{'ITG'}
% [js_elecNums,js_error,~,distToHippHead] = getNeuroLocLeads(subj,'ITG',[]);


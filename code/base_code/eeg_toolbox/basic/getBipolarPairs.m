function [allpairs, allTagNames, allGrpNames, allElecType] = getBipolarPairs(subj,source_subj)
% [allPairs, allTagNames, allGrpNames, allElecType] = ...
%      GETBIPOLARPAIRS(subj, source_subj)

MOUNT_DIR = '';
% directories
subjDir = fullfile(MOUNT_DIR,'/data/eeg',subj);
docDir  = fullfile(subjDir,'docs');
talDir  = fullfile(subjDir,'tal');

% file names
voxMotherFileName = 'VOX_coords_mother.txt';
jackFileName      = 'jacksheet.txt';
jackFileName_2    = 'jack_sheet.txt';
elecFileName      = 'electrodes.m';
leadsFileName     = 'leads.txt';


% get jacksheet
[elecNum elecNam elecNam_stripped elecNam_stripped_num] = ...
    getTheContentsOFTheJackFile_local(docDir,jackFileName,jackFileName_2);
if isempty(elecNum);error(' no jacksheet');end

% get the vox coords mother file
if exist('source_subj','var')&&~isempty(source_subj)
    talDir_source  = fullfile(MOUNT_DIR,'/data/eeg',source_subj,'tal');
    voxMomStruct = getTheContentsOFVoxMom_local(talDir_source,voxMotherFileName);
else
    voxMomStruct = getTheContentsOFVoxMom_local(talDir,voxMotherFileName);
end
if isempty(voxMomStruct); allpairs = nan;return;end

% get the electrodes.m file
gridsAndStrips = getTheContentsOFTheElecs_local(docDir,elecFileName);
allNeuralElecs = [];
for k=1:size(gridsAndStrips,1)
    allNeuralElecs = cat(1,allNeuralElecs,[gridsAndStrips(k,1):gridsAndStrips(k,2)]');
end

%check if it has in leads.txt
fid = fopen(fullfile(talDir,leadsFileName),'r');
c = textscan(fid,'%d');
leads = c{1};
clear c
fclose(fid)

% loop through the electrode groups
% make the bipolar montage for each one
unElecNam_stripped = unique(elecNam_stripped);
fprintf('\n')
fprintf('Looping all the electrode groups:\n')
allpairs    = [];
allTagNames = {};
allGrpNames = {};
allElecType = {};
for k=1:length(unElecNam_stripped)
    thisGrp    = unElecNam_stripped{k};
    fprintf('% 8.8s',thisGrp)
    thisGrpInd = find(strcmp(elecNam_stripped,thisGrp));
    thisGrpNam = elecNam(thisGrpInd);
    thisGrpNam_num = elecNam_stripped_num(thisGrpInd);
    thisGrpEls = elecNum(thisGrpInd);
    if sum(ismember(thisGrpEls,allNeuralElecs))==0
        fprintf(' not a neural channel...SKIPPING\n')
        continue
    end
    
    % get the geometry for this group
    numElsThisGroup = length(thisGrpNam);
    allGeometry     = [];
    throwTheseOut   = false(1,numElsThisGroup);
    all_sg_or_d     = {};
    %get geometry
    for k1=1:numElsThisGroup
        thisElecIndinVoxMom = find(strcmp({voxMomStruct.name},thisGrpNam{k1}));
        if isempty(thisElecIndinVoxMom)
            if ~ismember(thisGrpEls(k1),leads)
                % if it is not found in vox mom, but it is removed, then skip
                % it
                throwTheseOut(k1) = true;
                sprintf('Throwing out %s: found in jacksheet, not in vox mom',...
                    thisGrpNam{k1})
            else
                % if it is in leads.txt, but it is not found in vox mother,
                % then you want to throw an error because that is not what I
                % expect.
                error(sprintf('Found %s in jacksheet but not in vox mom',...
                    thisGrpNam{k1}))
            end
        elseif ~ismember(thisGrpEls(k1),leads)
             throwTheseOut(k1) = true;
                sprintf('Throwing out %s: not found in leads',...
                    thisGrpNam{k1})
        end
        allGeometry = cat(1,allGeometry, voxMomStruct(thisElecIndinVoxMom).geom);
        all_sg_or_d = cat(1,all_sg_or_d, voxMomStruct(thisElecIndinVoxMom).type);
    end
    
    unAllGeometry = unique(allGeometry,'rows');
    if size(unAllGeometry,1)>1
        fprintf('WARNING:::multiple geometries found for %s',thisGrp)
        %error('error HERE'); no error because this can occur if a grid is cut
        %into smaller subgrids, see below
    end
    unAll_sg_or_d = unique(all_sg_or_d);
    if size(unAll_sg_or_d,1)>1
        fprintf('this thing does not have a consistent grid or strip %s',thisGrp)
        error('error HERE')
    end
    
    %update these to result in a mismatch in geometry
    thisGrpInd(throwTheseOut) = [];
    thisGrpNam(throwTheseOut) = [];
    thisGrpEls(throwTheseOut) = [];
    
    % check to see the geometry makes sense
    numElsThisGroup = length(thisGrpNam);
    
    if numElsThisGroup==0
        % I REALLY HOPE THIS DOESN'T BREAK ANYTHING. -IMP 8/1/14
        continue
    end
    if numElsThisGroup~=prod(unAllGeometry)
        switch lower(upper(unAll_sg_or_d{1}))
            case {'s','d'}
                if length(thisGrpEls)<=1
                    continue
                end
                allpairs_tmp = [thisGrpEls(1:end-1) thisGrpEls(2:end)];
                
            case {'g'}
                if size(unAllGeometry,1) ==1 && numElsThisGroup< prod(unAllGeometry)
                    %This is to adjust for when contacts from the middle of the grid are not being recorded from
                    startVal                    = min(thisGrpEls);
                    
                    %This is if electrodes have not been skipped in the recording montage
                    if sum(diff(thisGrpEls)>1)==0 && sum(diff(thisGrpNam_num)>1)==0
                        
                        %compute allpairs_tmp so it ends in the right place.
                        %However, the pairs that occur before the break will have
                        %lower values than they should
                        allpairs_tmp  = make_a_bipolar_motage(startVal-(prod(unAllGeometry)-numElsThisGroup),...
                            unAllGeometry(1),...
                            unAllGeometry(2));
                        
                        %figure out which pairs are good
                        tmp = allpairs_tmp-min(min(allpairs_tmp))+1;
                        indsToKeep_tmp              = ismember(tmp,thisGrpNam_num);
                        
                        %loop through the inappropriate allpairs_tmp values
                        %and adjust appropriately
                        %find the break
                        diff_enam_num = diff(thisGrpNam_num);
                        breakLoc = find(diff_enam_num>1);
                        for bl = 1:length(breakLoc)
                            idx2adjust = tmp<=thisGrpNam_num(breakLoc(bl));   %all the elec nums that are lower than the break
                            allpairs_tmp(idx2adjust) = allpairs_tmp(idx2adjust) +  (diff_enam_num(breakLoc(bl))-1); %increase it by the number of elecs skipped
                        end
                    else  %This is when the missing leads are skipped in the recording montage
                        
                        if sum(diff(thisGrpEls)>1)>0 %skips are accounted for in jacksheet
                            allpairs_tmp                = make_a_bipolar_motage(startVal,...
                            unAllGeometry(1),...
                            unAllGeometry(2));
                            indsToKeep_tmp              = ismember(allpairs_tmp,thisGrpEls);
                        elseif  sum(diff(thisGrpNam_num)>1)>0   %skips not accounted for in jacksheet, only in elec names 
                            %make  pairs with elecGrpNam_num
                            thisStart =  min(thisGrpNam_num);
                            allpairs_namNum                = make_a_bipolar_motage(thisStart,...
                            unAllGeometry(1),...
                            unAllGeometry(2));
                            indsToKeep_tmp              = ismember(allpairs_namNum,thisGrpNam_num); % a bit messy because this gets re-done later
                            indsToKeep                  = sum(indsToKeep_tmp,2)==2;
                            allpairs_namNum(~indsToKeep,:) = [];
                            % populate allpairs_tmp based on namNum_pairs
                            for jj = 1:size(allpairs_namNum, 1)
                               allpairs_tmp(jj,1) = thisGrpEls(allpairs_namNum(jj,1)==thisGrpNam_num);
                               allpairs_tmp(jj,2) = thisGrpEls(allpairs_namNum(jj,2)==thisGrpNam_num);
                            end
                            indsToKeep_tmp = ones(size(allpairs_tmp,1),2);
                        end
                            
                    end
                    indsToKeep                  = sum(indsToKeep_tmp,2)==2;
                    allpairs_tmp(~indsToKeep,:) = [];
                else
                    %This is to adjust for when a grid has been cut into smaller subgrids or strips
                    %to do: account for when grid is cut into subgrids of varying
                    %sizes
                    allpairs_tmp = [];
                    for i = 1:size(unAllGeometry,1)
                        retain = ismember(allGeometry(:,1),unAllGeometry(i,1)) &...
                            ismember(allGeometry(:,2),unAllGeometry(i,2));
                        startVals = min(thisGrpEls(retain));
                        
                        thesePairs    = make_a_bipolar_motage(startVals,...
                            unAllGeometry(i,1),...
                            unAllGeometry(i,2));
                        
                        % Throw out pairs in which one electrode doesn't exist in
                        % this
                        thesePairs(~all(ismember(thesePairs, thisGrpEls),2),:) = [];
                        allpairs_tmp = [allpairs_tmp;thesePairs];
                    end
                    
                    %           startVals = min(thisGrpEls):prod(unAllGeometry):max(thisGrpEls);
                    %           allpairs_tmp = [];
                    %           for sv = 1:length(startVals)
                    %              thesePairs    = make_a_bipolar_motage(startVals(sv),...
                    %                               unAllGeometry(1),...
                    %                               unAllGeometry(2));
                    %               allpairs_tmp = [allpairs_tmp;thesePairs];
                    %           end
                end
            otherwise
                error('this is not a strip, grid, or a depth')
        end
        
    else
        startVal     = min(thisGrpEls);
        allpairs_tmp = make_a_bipolar_motage(startVal,unAllGeometry(1), ...
            unAllGeometry(2));
    end
    
    %convert
    [allnames_tmp] = pairs2names_local(allpairs_tmp,thisGrpEls,thisGrpNam);
    
    % update larger vars
    allpairs    = cat(1,allpairs,   allpairs_tmp);
    allTagNames = cat(1,allTagNames, allnames_tmp);
    allGrpNames = cat(1,allGrpNames,repmat({thisGrp},size(allpairs_tmp,1),1));
    allElecType = cat(1,allElecType,repmat(unAll_sg_or_d, ...
        size(allpairs_tmp,1),1));
    fprintf('\n')
end

save(fullfile('/data/eeg',subj,'/tal/bpPairs.mat'),'allpairs','allTagNames','allGrpNames','allElecType');


function out = getTheContentsOFVoxMom_local(d,f);
fil = fullfile(d,f);
fid = fopen(fil,'r');
if fid==-1
    out = [];
    fprintf('%s does not exist...EXITING\n\n',f)
    return
end
X = textscan(fid,'%s%d%d%d%s%s','delimiter','\t');
fclose(fid);
gridNames    = X{1};
elecXYZ      = [X{2} X{3} X{4}];
stripOrDepth = X{5};
geometry_tmp = X{6};
out          = [];
for k=1:size(geometry_tmp,1)
    out(k).name = gridNames{k};
    out(k).voxP = elecXYZ(k,:);
    out(k).type = stripOrDepth{k};
    out(k).geom = sscanf(geometry_tmp{k},'%d%d')';
end


function [num nam nam_stripped nam_strimmed_num] = getTheContentsOFTheJackFile_local(dDir,j1,j2);
jFil_1 = fullfile(dDir,j1);
jFil_2 = fullfile(dDir,j2);
fid_1  = fopen(jFil_1,'r');
fid_2  = fopen(jFil_2,'r');

if fid_1~=-1 & fid_2==-1
    %fprintf('Found jack sheet in %s\n\n',jFil_1)
    fid = fid_1;
elseif fid_1==-1 & fid_2~=-1
    %fprintf('Found jack sheet in %s\n\n',jFil_2)
    fid = fid_2;
elseif fid_1==-1 & fid_2==-1
    fprintf('Found no jack sheet\n\n')
    num = [];
    nam = [];
    nam_stripped=[];
    nam_strimmed_num=[];
    return
else
    error('CHECK THIS OUT....I HAVE NO IDEA WHAT THIS CONDITION MEANS')
end

X   = textscan(fid,'%f%s');
num = X{1};
nam = X{2};
clear X

Nelec       = size(nam,1);
elecJackTag = cell(Nelec,1);
elecJackNum = nan(Nelec,1);
for k=1:Nelec
    thisElecNam     = nam{k};
    elecJackTag{k}  = thisElecNam(regexp(thisElecNam,'\D'));
    elecJackNum_tmp = str2double(thisElecNam(regexp(thisElecNam,'\d')));
    elecJackNum(k)=elecJackNum_tmp;
end
nam_stripped = elecJackTag;
nam_strimmed_num = elecJackNum;

function out = getTheContentsOFTheElecs_local(d,f);
elecFile = fullfile(d,f);
if ~exist(elecFile,'file')
    out = [];
    fprintf('%s does not exist...EXITING\n\n',f)
    return
end
run(elecFile);
if exist('r','var')
    out = r;
elseif exist('grids','var')
    out = grids;
else
    out = [];
    fprintf('The contects of electrodes.m not right...EXITING\n\n')
    return
end

function  allpairs = make_a_bipolar_motage(startval,numElecs_dim1,numElecs_dim2)

% this is all Aswin G. Ramayya.  This occurred on June 14, 2013.
% Ashwin is completelty responsible for this code and miustakes
% therein.  John F. Burke is not reponsible for any errors,
% mistakes, or general bad things that come out of this code. But
% John Burke checked thouroughly and has placed his personal
% stamp of approval, thus accepting all responsibility (but no
% credit). Ashwin would like to add, spontaneously, that he is a
% tool. John agrees.
allpairs = [];
for j = 1:numElecs_dim2
    for i = 1:numElecs_dim1
        if i < numElecs_dim1
            allpairs(end+1,:) = [i-1 i] + startval + numElecs_dim1*((j-1));
        end
        if j < numElecs_dim2
            allpairs(end+1,:) = [i-1 i+numElecs_dim1-1] + startval + numElecs_dim1*((j-1));
        end
    end
end

function [allnames_tmp] = pairs2names_local(allpairs_tmp,thisGrpEls,thisGrpNam);
allnames_tmp = {};
for i = 1:size(allpairs_tmp,1)
    for j = 1:size(allpairs_tmp,2)
        allnames_tmp(i,j) = thisGrpNam(allpairs_tmp(i,j) == thisGrpEls);
    end
end
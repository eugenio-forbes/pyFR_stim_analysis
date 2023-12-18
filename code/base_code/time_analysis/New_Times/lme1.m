clearvars;
close all;
clc;
warning off;

load('hippoStruct.mat')
sesscodes = [hippoStruct(:).sesscode];
sesscodes = unique(sesscodes([hippoStruct(:).hasES]&[hippoStruct(:).hasNS]));

title = 'TimexCon_ESNS_0-2_AGB_slope'
formula = 'pow ~ time*condition + (1+condition|sub:elec) + (1|sub:trial)';
type = 'encoding';
norm = 'nostim';
variable = 'normPow';
conditions = {'ES','NS'};
times = [0,2];
ROIs = {'AGL','AGR'};
rootDir = '/project/TIBIR/Lega_lab/s427026/AGPCC';
tableDir = fullfile(rootDir,'Tables_WM_encoding',norm,variable);
saveDir =  fullfile(rootDir,'LMEM_time_new_23',type,norm,variable,title);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

AR_c = zeros(49,9);
AL_c = zeros(49,9);
PR_c = zeros(49,9);
PL_c = zeros(49,9);

for freq = 1:49
    
    file = fullfile(tableDir,num2str(freq),'t.mat');
    
    load(file)
    
    t = t(ismember(t.stimReg,ROIs),:);
    t = t(ismember(t.sess,sesscodes),:);
    t = t(ismember(t.time,times),:);
    t = t(ismember(t.condition,conditions),:);
    
    Rt = t(contains(t.stimReg,'R'),:);
    Lt = t(contains(t.stimReg,'L'),:);
    clear t
    L_AR = Lt(ismember(Lt.hipReg,'AR'),:);
    L_AL = Lt(ismember(Lt.hipReg,'AL'),:);
    L_PR = Lt(ismember(Lt.hipReg,'PR'),:);
    L_PL = Lt(ismember(Lt.hipReg,'PL'),:);
    clear Lt
    L_AR.hipReg = strrep(L_AR.hipReg,'R','L');
    L_AL.hipReg = strrep(L_AL.hipReg,'L','R');
    L_PR.hipReg = strrep(L_PR.hipReg,'R','L');
    L_PL.hipReg = strrep(L_PL.hipReg,'L','R');
    Lt = [L_AR;L_AL;L_PR;L_PL];
    clear L_AR L_AL L_PR L_PL
    t = [Rt;Lt];
    clear Rt Lt
    
    t.sub = categorical(t.sub);
    t.sess = categorical(t.sess);
    t.elec = categorical(t.elec);
    t.condition = contains(t.condition,'ES');
    t.time = t.time == 2;
    
    AR = t(ismember(t.hipReg,'AR'),:);
    AL = t(ismember(t.hipReg,'AL'),:);
    PR = t(ismember(t.hipReg,'PR'),:);
    PL = t(ismember(t.hipReg,'PL'),:);
    
    AR_lme = fitlme(AR,formula);
    AL_lme = fitlme(AL,formula);
    PR_lme = fitlme(PR,formula);
    PL_lme = fitlme(PL,formula);
    
    AR_c(freq,1:3) = [AR_lme.Coefficients.tStat(2:4)];
    AL_c(freq,1:3) = [AL_lme.Coefficients.tStat(2:4)];
    PR_c(freq,1:3) = [PR_lme.Coefficients.tStat(2:4)];
    PL_c(freq,1:3) = [PL_lme.Coefficients.tStat(2:4)];
    AR_c(freq,4:6) = [AR_lme.Coefficients.pValue(2:4)];
    AL_c(freq,4:6) = [AL_lme.Coefficients.pValue(2:4)];
    PR_c(freq,4:6) = [PR_lme.Coefficients.pValue(2:4)];
    PL_c(freq,4:6) = [PL_lme.Coefficients.pValue(2:4)];
    AR_c(freq,7:9) = [AR_lme.Coefficients.Estimate(2:4)];
    AL_c(freq,7:9) = [AL_lme.Coefficients.Estimate(2:4)];
    PR_c(freq,7:9) = [PR_lme.Coefficients.Estimate(2:4)];
    PL_c(freq,7:9) = [PL_lme.Coefficients.Estimate(2:4)];
    clear t Enc Rec AR AL PR PL AR_lme AL_lme PR_lme PL_lme
end
        
save(fullfile(saveDir,'data.mat'),'AL_c','AR_c','PL_c','PR_c')
save(fullfile(saveDir,'formula.mat'),'formula')
save(fullfile(saveDir,'ROIs.mat'),'ROIs')
save(fullfile(saveDir,'conditions.mat'),'conditions')
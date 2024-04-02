function subj_match_EDB = getSubjDominHem(hem,disp)

if ~exist('disp','var')||isempty(disp)
  disp=true;
end

% get the data from the look up table(LUT)
LUT = 'basicECoGSubjInfo.txt';
fid = fopen(LUT);
INFO = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s');
subj_LUT      = INFO{:,1}; 
hand_LUT      = upper(INFO{:,2}); 
gender_LUT    = upper(INFO{:,3}); 
age_LUT       = upper(INFO{:,4}); 
iq_LUT        = upper(INFO{:,5}); 
ethnicity_LUT = upper(INFO{:,6}); 
language_LUT  = upper(INFO{:,7});
langDom_LUT   = upper(INFO{:,8}); 

% get the subjects from the event database(EDB)
newSubjStr = 'preliminary data; not Talairached; not properly re-referenced';
badLeadStr = 'no bad leads info';

subj_EDB_all  = get_subs('pyFR');
subj_EXPINFO  = get_sub_expInfo('pyFR',subj_EDB_all);
subj_good     = [subj_EXPINFO.hadProblems]==0;
subj_newer    = strcmp({subj_EXPINFO.description},badLeadStr);
subj_notready = strcmp({subj_EXPINFO.description},newSubjStr);

subj_EDB = subj_EDB_all(subj_good|subj_newer);

% see if you have all the subjects info
missing_subj = subj_EDB(~ismember(subj_EDB,subj_LUT));
if ~isempty(missing_subj)
  if disp
    fprintf('\nWARNING: Missing info from the following subjects:\n')
    for k=1:length(missing_subj)
      fprintf('   %s\n',missing_subj{k})
    end
    fprintf('\ncontinuing...\n\n',missing_subj{k})
  end
end

% get all subjects in the LUT that match this request
hem = upper(hem);
switch hem
 case 'RIGHT'
  ind = find(strcmp(hand_LUT,'LEFT')|strcmp(langDom_LUT,'RIGHT'));
 case 'LEFT'
  ind = find(strcmp(hand_LUT,'RIGHT')|strcmp(langDom_LUT,'LEFT'));
 case 'NOINFO'
  ind = find(strcmp(hand_LUT,'-999')&strcmp(langDom_LUT,'-999'));
 case 'ALL'
  subj_match_EDB = subj_EDB;
  return
 otherwise
  error('hem must be RIGHT, LEFT or NoInfo')
end
subj_match_LUT = subj_LUT(ind);

% now only output those subjects that are also in the EDB
subj_match_EDB = subj_EDB(ismember(subj_EDB,subj_match_LUT));
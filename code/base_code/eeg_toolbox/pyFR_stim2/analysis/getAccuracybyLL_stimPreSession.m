function getAccuracybyLL_stimPreSession(ev,doConds)
%
% FUNCTION:
%   
% DESCRIPTION:
%
% INPUT:
%
% OUTPUT:
%
% NOTES:
%   (1) written by jfburke [''] (john.fred.burke@gmail.com)
%
%



% start with the words
ev_all = ev;
ev     = ev(strcmp({ev.type},'WORD'));

% ste up vectors
allLL     = unique([ev.serialpos]);
accByLL   = cell(1,length(allLL));
accByLL_h = cell(1,length(allLL));
accByLL_l = cell(1,length(allLL));

unSess = unique([ev.session]);
for s=1:length(unSess)
  thisSess = unSess(s);
  ev_sess  = ev([ev.session]==thisSess);
  
  unList = unique([ev_sess.list]);
  for s=1:length(unList)
    thisList    = unList(s);    
    ev_sessList =  ev_sess([ev_sess.list]==thisList);
    
    thisLL          = max([ev_sessList.serialpos]);
    thisAcc         = sum([ev_sessList.recalled]==1)./length(ev_sessList);
    LLind           = find(allLL==thisLL);
    accByLL{LLind}  = cat(1,accByLL{LLind},thisAcc);
    
    if strcmp(upper(ev_sessList(1).imagery),'HIIMAG')
      accByLL_h{LLind}  = cat(1,accByLL_h{LLind},thisAcc);
    else
      accByLL_l{LLind}  = cat(1,accByLL_l{LLind},thisAcc);
    end
        
  end  
end

% remove empty
emptyInd            = cellfun(@isempty,accByLL);
accByLL(emptyInd)   = [];
accByLL_h(emptyInd) = [];
accByLL_l(emptyInd) = [];
allLL(emptyInd)     = [];

% plot the curve
if ~doConds
  ave = nan(size(allLL));
  sem = nan(size(allLL));
  for k=1:length(allLL)
    ave(k) = mean(accByLL{k});
    sem(k) = std(accByLL{k})./sqrt(length(accByLL{k}));  
  end
  errorbar(allLL,ave,sem,'LineWidth',3,'Color','k')
else
  hold on
  ave = nan(size(allLL));
  sem = nan(size(allLL));
  for k=1:length(allLL)
    ave(k) = mean(accByLL_h{k});
    sem(k) = std(accByLL_h{k})./sqrt(length(accByLL_h{k}));  
  end
  errorbar(allLL,ave,sem,'LineWidth',3,'Color','b')
  
  
  ave = nan(size(allLL));
  sem = nan(size(allLL));
  for k=1:length(allLL)
    ave(k) = mean(accByLL_l{k});
    sem(k) = std(accByLL_l{k})./sqrt(length(accByLL_l{k}));  
  end
  errorbar(allLL,ave,sem,'LineWidth',3,'Color','r')
  
  legend('HI Imagery','LO Imagery')
  hold off
end
    
formatPic(gca)
grid on
box off
xlabel('list length')
ylabel('Probability of recall')
title('Baseline Memory Task')


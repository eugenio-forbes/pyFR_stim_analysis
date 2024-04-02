subj = {'CC002' 'CC007' 'CC009' 'CC013' 'CC014' 'CC015' 'CC016' 'CC020' 'CC024' 'UT001' 'UT005' 'UT008' 'UT010' 'UT012'};

windWidth = 50;
numWinds = 900/50;

for n = 1:length(subj)
    clear zMatNonCollapsed
    clear zMa
    subjDir= sprintf('/Users/brad/Documents/HippSync/%s/psync250',subj{n});
    load(strcat(subjDir,'/APHippSynch.mat'));
    if n == 1
        %         a = mean(mean(zMatrixNon));
        %         b = mean(mean(zMatrixRec));
        zMatNon = reshape(zMatrixNon,[pairCounter,45,900]);
        zMatRec = reshape(zMatrixRec,[pairCounter,45,900]);
        
        for wind_ind = 1:numWinds
            %Find the index of the first sample for this window
            windStart = 1 + (wind_ind-1)*windWidth;
            %Find the index of the last sample for this window
            windEnd = windStart + windWidth - 1;
            zMatNonCollapsed(:,:,wind_ind) = nanmean(zMatNon(:,:,windStart:windEnd),3);
            zMatRecCollapsed(:,:,wind_ind) = nanmean(zMatRec(:,:,windStart:windEnd),3);
        end
        
        zMatrixNonAllCollapsed = zMatNonCollapsed;
        zMatrixRecAllCollapsed = zMatRecCollapsed;
        totalElec = pairCounter;
    else
        %         c = mean(mean(zMatrixNon));
        %         d = mean(mean(zMatrixRec));
        %         a = a + c;
        %         b = b + d;
        
        zMatNon = reshape(zMatrixNon,[pairCounter,45,900]);
        zMatRec = reshape(zMatrixRec,[pairCounter,45,900]);
        for wind_ind = 1:numWinds
            %Find the index of the first sample for this window
            windStart = 1 + (wind_ind-1)*windWidth;
            %Find the index of the last sample for this window
            windEnd = windStart + windWidth - 1;
            zMatNonCollapsed(:,:,wind_ind) = nanmean(zMatNon(:,:,windStart:windEnd),3);
            zMatRecCollapsed(:,:,wind_ind) = nanmean(zMatRec(:,:,windStart:windEnd),3);
        end
        zMatrixNonAllCollapsed = vertcat(zMatrixNonAllCollapsed,zMatNonCollapsed);
        zMatrixRecAllCollapsed = vertcat(zMatrixRecAllCollapsed,zMatRecCollapsed);
        
        totalElec = totalElec + pairCounter;
    end
end

%Avg stuff
a = mean(zMatrixNonAllCollapsed);
b = mean(zMatrixRecAllCollapsed);

NonAvg = squeeze(a);
RecAvg= squeeze(b);

diff = RecAvg-NonAvg;

%ttest stuff
[~,p,~,stat] = ttest2(zMatrixRecAllCollapsed,zMatrixNonAllCollapsed);
           
t = squeeze(stat.tstat);
p = squeeze(p);

%Background plotting settings
Rs = 1000;

nyquist_freq = floor(Rs/2);

freqs = eeganalparams('freqs');

valid_freqs_ind = find(freqs<nyquist_freq);

valid_freqs = freqs(valid_freqs_ind);
freqs = valid_freqs;
freqs = freqs(1:45);

encodeOffset = 0;
encodeDur = 1800;
timeVect = (encodeOffset+1):2:(encodeOffset+encodeDur);
timeVect = timeVect(1:50:length(timeVect));

%wait before plotting
keyboard
           
           figure(1)
           subplot(2,2,1) 
           contourf(timeVect,freqs,NonAvg,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('Electro Pair Avg NonRec');
           title(title_string);
           caxis([0 3])
           colorbar

           subplot(2,2,2)           
           contourf(timeVect,freqs,RecAvg,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('Electro Pair Avg Rec');
           title(title_string);
           caxis([0 3])
           colorbar
           
           figure(2)
           contourf(timeVect,freqs,diff,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('Avg Rec-Non');
           title(title_string);
           caxis([-0.6 0.7])
           colorbar
           
           figure(3)
           subplot(2,2,1) 
           contourf(timeVect,freqs,t,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('ttest');
           title(title_string);
           caxis([-4 4])
           colorbar

           subplot(2,2,2)           
           contourf(timeVect,freqs,p,40,'edgecolor','none')
           set(gca,'yscale','log','ytick',[2.^(1:8)]);
           title_string = sprintf('p-values');
           title(title_string);
           caxis([0 1])
           colorbar
           
           
           
%         c = zMatrixNon;
%         d = zMatrixRec;
%         a = a + c;
%         b = b + d;
%        
%            
%         zMatrixNonVec = reshape(zMatrixNon,[pairCounter,45,900]);
%         zMatrixRecVec = reshape(zMatrixRec,[pairCounter,45,900]);
%            
%         zMatrixNonAll = vertcat(zMatrixNonAll,zMatrixNonVec);
%         zMatrixRecAll = vertcat(zMatrixRecAll,zMatrixRecVec);

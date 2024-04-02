function update_tal_isGood(subject)
%UPDATE_TAL_ISGOOD Update Talairach database with new good_leads.txt info.
%   If a patient has already been added to the Talairach database and
%   you've since gotten new bad leads info on that patient, use this
%   function to update the Talairach database based on the patient's
%   good_leads.txt file. A backup copy of the file containing Talairach
%   information for all patients will be created each time changes are made
%   to it (in case the changes need to be reverted for any reason).
%   
%   example:
%   update_tal_isGood('TJ037')
%
%   Written by Ryan Williams (March 2012).
%   Last updated on 2012-03-27.
    
    %load the current Talairach database (file containing all patients)
    talDatabaseFilename = '/data/eeg/tal/allTalLocs_GM.mat';
    load(talDatabaseFilename)
    
    %load good_leads.txt for the specified patient
    goodLeadsFilename = ['/data/eeg/', subject, '/tal/good_leads.txt'];
    goodLeads = load(goodLeadsFilename);
    leadsIndicesTal = find(strcmp(subject, {events.subject}));
    
    %note the year (for making a backup of the Talairach database)
    currentClock = clock;
    backupFilename = [talDatabaseFilename, '.', num2str(currentClock(1)), '-'];
    
    %note the month
    if currentClock(2) < 10
        backupFilename = [backupFilename, '0'];
    end
    backupFilename = [backupFilename, num2str(currentClock(2)), '-'];
    
    %note the day of the month
    if currentClock(3) < 10
        backupFilename = [backupFilename, '0'];
    end
    backupFilename = [backupFilename, num2str(currentClock(3)), '_'];
    
    %note the hour
    if currentClock(4) < 10
        backupFilename = [backupFilename, '0'];
    end
    backupFilename = [backupFilename, num2str(currentClock(4)), '-'];
    
    %note the minute
    if currentClock(5) < 10
        backupFilename = [backupFilename, '0'];
    end
    backupFilename = [backupFilename, num2str(currentClock(5)), '-'];
    
    %note the second
    if currentClock(6) < 10
        backupFilename = [backupFilename, '0'];
    end
    backupFilename = [backupFilename, num2str(floor(currentClock(6)))];
    
    %make a backup of the current Talairach database (file containing all patients)
    copyfile(talDatabaseFilename, backupFilename);
   
    %update "isGood" field in the Talairach database (file containing all patients) where necessary
    numChangesFullDatabase = 0;
    for curTalIndexNum = 1 : length(leadsIndicesTal)
        isGood = 0;
        for curLeadsNum = 1 : length(goodLeads)
            if events(leadsIndicesTal(curTalIndexNum)).channel == goodLeads(curLeadsNum)
                isGood = 1;
            end
        end
        if events(leadsIndicesTal(curTalIndexNum)).isGood ~= isGood
            events(leadsIndicesTal(curTalIndexNum)).isGood = isGood;
            numChangesFullDatabase = numChangesFullDatabase + 1;
        end
    end
    
    %save any changes to the Talairach database (file containing all patients); alternatively, if no changes were made, delete the backup file
    if numChangesFullDatabase ~= 0
        save(talDatabaseFilename, 'events')
    else
        delete(backupFilename);
    end
    
    %load the current Talairach database (file containing a single patient)
    talDatabaseThisPatientFilename = ['/data/eeg/tal/allTalLocs_GM_', subject, '.mat'];
    load(talDatabaseThisPatientFilename)
    
    %update "isGood" field in the Talairach database (file containing a single patient) where necessary
    numChangesThisPatientDatabase = 0;
    for curTalIndexNum = 1 : length(events)
        isGood = 0;
        for curLeadsNum = 1 : length(goodLeads)
            if events(curTalIndexNum).channel == goodLeads(curLeadsNum)
                isGood = 1;
            end
        end
        if events(curTalIndexNum).isGood ~= isGood
            events(curTalIndexNum).isGood = isGood;
            numChangesThisPatientDatabase = numChangesThisPatientDatabase + 1;
        end
    end
    
    %save any changes to the Talairach database (file containing a single patient)
    if numChangesThisPatientDatabase ~= 0
        save(talDatabaseThisPatientFilename, 'events')
    end
    
    if numChangesFullDatabase == 0 && numChangesThisPatientDatabase == 0
        fprintf('\nWarning: Nothing was changed in %s nor %s.\nHas %s changed since the last time the Talairach database was updated?\n(If concerned, check backup versions of %s.)\n\n', talDatabaseFilename, talDatabaseThisPatientFilename, goodLeadsFilename, talDatabaseFilename);
    else
        %warn if no "isGood" values changed (file containing all patients), otherwise indicate how many changes were made
        if numChangesFullDatabase == 0
            fprintf('\nWarning: Nothing was changed in %s.\nHas %s changed since the last time the Talairach database was updated?\n(If concerned, check backup versions of %s.)\n\n', talDatabaseFilename, goodLeadsFilename, talDatabaseFilename);
        elseif numChangesFullDatabase == 1
            fprintf('\nUpdated 1 electrode in %s.\n\n', talDatabaseFilename);
        else
            fprintf('\nUpdated %d electrodes in %s.\n\n', numChangesFullDatabase, talDatabaseFilename);
        end

        %warn if no "isGood" values changed (file containing a single patient), otherwise indicate how many changes were made
        if numChangesThisPatientDatabase == 0
            fprintf('Warning: Nothing was changed in %s.\nHas %s changed since the last time the Talairach database was updated?\n(If concerned, check backup versions of %s.)\n\n', talDatabaseThisPatientFilename, goodLeadsFilename, talDatabaseFilename);
        elseif numChangesThisPatientDatabase == 1
            fprintf('Updated 1 electrode in %s.\n\n', talDatabaseThisPatientFilename);
        else
            fprintf('Updated %d electrodes in %s.\n\n', numChangesThisPatientDatabase, talDatabaseThisPatientFilename);
        end
    end
end


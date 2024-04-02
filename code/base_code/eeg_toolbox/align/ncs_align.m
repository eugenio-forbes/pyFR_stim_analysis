function[events] = ncs_align(sync_file,beh_file,events_file,chan_file)
%ncs_align  Compute eeg offsets for the given event times using the
%           pulse times from the behavioral and neuralynx files.

% 2-24-10  JRM  Wrote it.

if ~exist('chan_file','var')
    chan_file = fullfile(fileparts(sync_file),'chan');
end

minPulses = 200; %minimum number of pulses required for alignment
threshMS = 100; %maximum allowable deviation 

%load in pulses
[~,pulse_str] = system(['cut -f 1 ' sync_file]);
pulse_samples = strread(pulse_str,'%u'); %sample numbers of received pulses

pulse_times = textread(beh_file,'%f\t%*d\t%*s\n'); %behavioral times pulses were sent (unix time, msec)

assert(length(pulse_samples) >= length(pulse_times),...
       'not all sent pulses were received!');

%load events struct
events = load(events_file); events = events.events;
beh_times = [events.mstime]; %times of behavioral events (unix time, msec)

%normalize pulse_times and beh_times using min([pulse_times beh_times])
first_time = min([pulse_times ; beh_times']);
pulse_times = pulse_times - first_time;
beh_times = beh_times - first_time;

errs = nan(4,length(pulse_samples)-minPulses);
bs = nan(size(errs,2),2);
for i = 1:length(errs)    
    [next_pulse_samples,next_pulse_times] = trim_to_length(pulse_samples,pulse_times,i);
    
    %compute regression line between pulse sample numbers and behavioral times
    bs(i,:) = regress(next_pulse_samples,...
                      [next_pulse_times ones(size(next_pulse_times))]);

    next_estimated_samples = bs(i,1)*next_pulse_times + bs(i,2);
    errs(1,i) = max(abs(next_estimated_samples - next_pulse_samples));
    errs(2,i) = median(abs(next_estimated_samples - next_pulse_samples));
    errs(3,i) = mean(abs(next_estimated_samples - next_pulse_samples));
    errs(4,i) = min(abs(next_estimated_samples - next_pulse_samples));
    subplot(4,1,1);
    plot(0:size(errs,2)-1,errs(1,:),'k','LineWidth',2); ylabel('Max'); title('Alignment Error');
    subplot(4,1,2);
    plot(0:size(errs,2)-1,errs(2,:),'k','LineWidth',2); ylabel('Median');
    subplot(4,1,3);
    plot(0:size(errs,2)-1,errs(3,:),'k','LineWidth',2); ylabel('Mean');
    subplot(4,1,4);
    plot(0:size(errs,2)-1,errs(4,:),'k','LineWidth',2); ylabel('Min'); xlabel('Number of pulses skipped');
    drawnow;
end

best_ind = find(errs(1,:) == min(errs(1,:)));
assert(length(best_ind) == 1,'multiple alignments were found!');
assert(min(errs(1,:)) <= threshMS,'good alignment not found!');
[pulse_sample_nums,pulse_times_sent] = trim_to_length(pulse_samples,pulse_times,best_ind);
eeg_offsets = round(interp1(pulse_times_sent,pulse_sample_nums,beh_times));
assert(length(events) == length(eeg_offsets),'not all events have been assigned an offset!');

for i = 1:length(events)
    if isnan(eeg_offsets(i))
        events(i).eegfile = '';
        events(i).eegoffset = 0;
    else
        events(i).eegfile = chan_file;
        events(i).eegoffset = eeg_offsets(i);
    end
end

try
    copyfile(events_file,[events_file,'.old']);
catch %#ok<CTCH>
    unix(['cp ',events_file,' ',events_file,'.old']);
end
save(events_file,'events');

function[a,b] = trim_to_length(a,b,i)
a = a(i:end);

if length(a) > length(b)
    a = a(1:length(b));
elseif length(b) > length(a)
    b = b(1:length(a));
end


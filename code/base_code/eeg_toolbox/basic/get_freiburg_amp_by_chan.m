function[amp] = get_freiburg_amp_by_chan(chan)

N_CHANS_PER_AMP = 8;

amp = floor((chan-1)/N_CHANS_PER_AMP) + 1;
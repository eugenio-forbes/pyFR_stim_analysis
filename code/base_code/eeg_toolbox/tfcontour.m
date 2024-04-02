function h = tfcontour(dat,freqs,durationMS,offsetMS,numContours,frange,trangeMS)
%TFCONTOUR - Wrapper to contourf to make a nice TF-plot.
%
% FUNCTION:
%   h = tfcontour(dat,freqs,durationMS,offsetMS,numContours,frange,trangeMS)
%
% INTPUT ARGS:
%
%
% OUTPUT ARGS:
%
%

if ~exist('numContours','var') | isempty(numContours)
  numContours = 25;
end
if ~exist('frange','var') | isempty(frange)
  frange = [min(freqs) max(freqs)];
end
if ~exist('trangeMS','var') | isempty(trangeMS)
  trangeMS = [offsetMS offsetMS+durationMS];
end

% set the ranges
% based on duration, offset, and samplerate
samplerate = fix(1000*size(dat,2)/durationMS);
offset = fix(offsetMS*samplerate/1000);
trange = fix(trangeMS*samplerate/1000);
trange = trange - offset + 1;
trange(2) = trange(2) - 1;
tbin = trange(1):trange(2);

% get the fbin
fbin = [min(find(freqs>=frange(1))):max(find(freqs<=frange(2)))];
freqs = freqs(fbin);

% resize the dat
dat = dat(fbin,tbin);

x = linspace(trangeMS(1),trangeMS(2),size(dat,2));
y = freqs;
fticks = 2.^[0:7];
fticks = fticks(fticks>=min(freqs) & fticks<=max(freqs));
h = contourf(x,y,double(dat),numContours,'linestyle','none');
set(gca,'yscale','log','ytick',fticks);
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
colorbar











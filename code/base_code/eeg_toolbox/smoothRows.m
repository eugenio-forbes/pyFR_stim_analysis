	
function y=smoothRows(x,windSize)
if ~exist('windSize','var')
 windSize=40; %at a 500-hz sampling rate, this is a 80 ms window (80 ms / .02s sampling period)
end

y=convn(x,ones(1,1,windSize)./windSize,'same');
function amp=normalizeBySessions(amp,events)
assert(size(amp,1)==length(events));

grp=[events.session]';
uGrp=unique(grp);
for i=1:length(uGrp)
  ind=grp==uGrp(i);
  m=amp(ind,:);
  me=mean(m(:));
  s=std(m(:));
  amp(ind,:)=(m-me)./s;
end

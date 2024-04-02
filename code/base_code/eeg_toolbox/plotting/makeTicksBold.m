% set the tics on all children
ticFontSize = 20;
hc = get(gcf,'Children');
for j=1:length(hc)
  if isprop(hc(j),'FontSize') & isprop(hc(j),'FontWeight')
    set(hc(j),'FontSize',ticFontSize,'FontWeight','Bold');
  end
end

function make_colored_surface(surf, annot)
%MAKE_COLORED_SURFACE(surf, annot) makes a colored surface based off of the
%   annotation information available in the annot file

[v, f] = read_surf(surf);
[~,label, colortable] = read_annotation(annot);
names = {};
for i = 1:size(colortable.table,1)
    this_label = colortable.table(i, 5);
    indices = (1:size(v));
    indices(label~=this_label)=[];
    f_label = f(any(ismember(f, indices),2),:);
    if size(f_label,1)>0
        hold all
        plotsurf_wrapper(v, f_label+1, colortable.table(i,1:3)./256,0);
        names = [names; colortable.struct_names{i}];
    end
end
legend(names)


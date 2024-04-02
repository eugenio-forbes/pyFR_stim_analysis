function plotWithLimits( handles )
%PLOTWITHLIMITS plots the 3d CT with limits on what is visible in the x, y,
%and z directions

limits = nan(3,2);
sliderHandles = [handles.xSlider1, handles.xSlider2; ...
            handles.ySlider1, handles.ySlider2;...
            handles.zSlider1, handles.zSlider2];
%loop over the slider handles, storing their percentages*600 to the limits.
%TODO: Change this so that it actually takes into account the max size of
%the pointcloud, rather than just assuming 600
for i=1:2
    for j=1:3
        sliderHandle = sliderHandles(j,i);
        val = get(sliderHandle,'value');
        min = get(sliderHandle,'min');
        max = get(sliderHandle,'max');
        percent = (val-min)./(max-min);
        limits(j,i)=percent.*600;
    end
end

%Plot with limits, keep the figure the same size
clickA3DPoint(handles.xyz', handles, limits);
xlim([0,500]);
ylim([0,500]);
zlim([0,600]);
h = get(handles.mainAxes,'children');
h_tags = get(h,'tag');

%put everything that isn't the fullBrain on top.
for t=1:length(h_tags)
    tag = h_tags{t};
    if ~strcmp(tag,'fullBrain')
        uistack(h(t),'top')
    end
end


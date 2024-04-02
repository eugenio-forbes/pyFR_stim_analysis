function[h]=plotSEM(x,y,SEM,kolor)
%This function plots "continuous errorbars."
%Inputs
%x                  x-values to plot
%y                  y-values to plot
%SEM                height of error bars (std. error of mean or 95 conf
                    %int)
%kolor              color of surface

%Output
%h                  handle of fill object

%Written by AGR ~10/2011. Updated 10/2012
if ~exist('kolor','var') || isempty(kolor)
    kolor = [0.01 0.01 0.99]; % default color is blue
end

h=fill([x,flipdim(x,2)],[y+SEM,flipdim(y-SEM,2)],kolor);
set(h,'edgecolor',kolor)

function [Avg, SEM] = StdErr(A)

%Function helps calculate the mean and standard error of mean of a
%distribution. Mainly used when you're trying to plot errorbars
%
%A = vector of the distribution you want the output on
%
%Last Updated: 08-31-2017
%Author: Jui-Jui Lin 


STD = nanstd(A,[],1);
SEM = STD ./ sqrt(size(A,1));
Avg = mean(A);
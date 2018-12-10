%% openNMS Matlab demonstration
% This script demonstrates the use of openNMS matlab function and shows
% how to display a 3D NMS Dataset nicely in matlab.

% Copyright 2016, Dr. Georg Wiora, NanoFocus AG

% Set input filename
inputFileName = '1-euro-star.nms';
% Set title of plot
plotTitle = 'Demonstration of NanoFocus openNMS';

% Run the plot function to read and display NMS-file
PlotNMS3D(inputFileName,plotTitle);

function hWaitFig=waitbarTimeRemaining2(waitbarProgress, waitBarTitle)
%function waitbarTimeRemaining(waitbarProgress)
%
% Functional purpose: To dispaly a modified waitbar with values of "ealpsed
%   time" and "remaining time"
%
% Input arguments:
%   hWaitbar- a handle to matlab waitbar (it is assumed that appropriate title is added to the waibar figure.
%
%   hTic- a handle to a tic enabled at the begining of the measured processes.
%
%   waitbarProgress- value between [0,1] describing procces done vs total process that is- current index/total indexes
%
% Output Arguments: None
%
% Issues & Comments:
%
% Author and Date:  Nikolay Skarbnik 14/04/2011
% Last update: 18/12/2012  waitbarProgress casted to double

if nargin<1
    waitbarProgress=0;
end

if nargin<2
    waitBarTitle='waitbarTimeRemaining';
    % [StackData, iStack]=dbstack(1);
    % waitBarTitle=StackData.name;
end

persistent hWaitbar;
persistent hTic;

if isempty(hWaitbar)
    % callingFunction=evalin('caller', 'mfilename');
    hWaitbar=waitbar(0, waitBarTitle, 'Name', waitBarTitle);
    set(hWaitbar, 'CloseRequestFcn',@waitbarProgressCloseFcn)
end

if isempty(hTic)
    hTic=tic;
end

if waitbarProgress>=1
    waitbarProgressCloseFcn;
%     close(hWaitbar);
%     clear('hWaitbar', 'hTic');
%     return;
end

if isempty(hTic) || isempty(hWaitbar)
    return;
end

tocVal=toc(hTic)*1e-5;
waitbarProgress=max( eps, double(waitbarProgress) ); % for some reason datestr fails to work with single type

waitbar(waitbarProgress, hWaitbar,...
    sprintf('Time  passed   %s  [H:Min:Sec.mSec].\nTime remaining %s  [H:Min:Sec.mSec].',...
    datestr(tocVal, 'HH:MM:SS.FFF'),...
    datestr(((1-waitbarProgress)/waitbarProgress*tocVal),'HH:MM:SS.FFF')));
if ~strcmpi(waitBarTitle, 'waitbarTimeRemaining') && ~strcmpi( waitBarTitle, get(hWaitbar, 'Name') )
    set(hWaitbar, 'Name', waitBarTitle);
end

    function waitbarProgressCloseFcn(~, ~) % src, evnt
        delete(hWaitbar);
        hWaitbar=[];
        hTic=[];
        %clear('hWaitbar', 'hTic');
        % close(hWaitbar);
    end
hWaitFig=hWaitbar;
end


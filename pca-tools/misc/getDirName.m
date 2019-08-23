function [fout, dout] = getDirName(din)
%% getDirName: function to parse folder name from current path
% This is simple really. In fact it's just a wrapper for the built-in fileparts
% function. I didn't know that function before I made this one. So yea.
%
% Usage:
%   [fout, dout] = getDirName(din)
%
% Input:
%   name_in:
%
% Output:
%   name_out:
%
%

% if isunix
%     p = regexp(name_in, '\/', 'split');
%     
% else
%     p = regexp(fin, '\\', 'split');
% end
% fout = p{end};

[dout, fout, ext] = fileparts(din);
fout              = [fout ext];

end
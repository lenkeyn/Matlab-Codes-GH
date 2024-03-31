%
% Use 'dir' as extension to filter for folders
%
function [files, fullPath] = findFileInSubfolders(Filepath, expression, ext)


if nargin < 3
    ext = [];
end

d = dir([Filepath, '\**\*']);

j = 1;
for i = 3:1:length(d)
    nMatch = 0;
    [~, ~, extention] = fileparts(d(i).name);
    if iscell(expression)
        for e = 1:1:numel(expression)
            if isequal(ext, 'dir')
                if d(i).isdir == true
                    nMatch = nMatch + numel(strfind(d(i).name, expression{e}));
                end
            elseif isequal(ext, []) || numel(strfind(extention, ext)) > 0
                nMatch = nMatch + numel(strfind(d(i).name, expression{e}));
            end
        end
    else
        if isequal(ext, 'dir')
            if d(i).isdir == true
                nMatch = nMatch + numel(strfind(d(i).name, expression));
            end
        elseif isequal(ext, []) || numel(strfind(extention, ext)) > 0
            nMatch = nMatch + numel(strfind(d(i).name, expression));
        end
    end

    if nMatch > 0
        files(j) = {d(i).name};
        fullPath(j) = {fullfile(d(i).folder,d(i).name)};
        j = j+1;
    end
end

%if length(files) == 1
%    files = files{1};
%    fullPath = fullPath{1};
%end

end
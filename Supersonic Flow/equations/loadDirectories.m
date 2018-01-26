function loadDirectories(varargin)

if ~isempty(varargin) && length(varargin)>1
    wdir = [pwd '/' varargin{2}];
else
    wdir = pwd;
end
folderNames = {'boundary', 'postProcess','solvers','visc','equations'};

cd ../
for i = 1:length(folderNames)
    addpath([pwd '/' folderNames{i}]);
end

cd(wdir)
if ~isempty(varargin)
    for i = 1:length(varargin{1})
        addpath([pwd '/' varargin{1}{i}]);
    end
end

end
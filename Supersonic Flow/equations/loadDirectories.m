function loadDirectories(eqnDir)

folderNames = {'boundary', 'postProcess','solvers','visc','equations'};
wdir = [pwd '/' eqnDir];

cd ../
for i = 1:length(folderNames)
    addpath([pwd '/' folderNames{i}]);
end

cd(wdir)

end
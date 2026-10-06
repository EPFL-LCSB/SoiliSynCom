function addSoilSynComDependencies(paths)
%ADDSOILSYNCOMDEPENDENCIES Add configured external toolboxes to MATLAB.

if nargin < 1
    paths = soilSynComPaths();
end

names = fieldnames(paths.external);
for index = 1:numel(names)
    location = paths.external.(names{index});
    if ~isempty(location) && isfolder(location)
        addpath(genpath(location));
    end
end
end

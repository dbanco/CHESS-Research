function outputs = loadOutputFile(inDir,inds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% === Build prefix ===
prefix = sprintf('output_j%d_%d_%d_', inds(1), inds(2), inds(3));

% === Get list of all PNG files ===
files = dir(fullfile(inDir, '*.mat'));

% === Find matching file ===
match_idx = [];
for k = 1:length(files)
    if startsWith(files(k).name, prefix)
        match_idx = [match_idx, k];
    end
end

% === Make sure there is exactly 1 match ===
if numel(match_idx) ~= 1
    error('Expected exactly one matching file. Found %d matches.', numel(match_idx));
end

% === Load the image ===
file_path = fullfile(inDir, files(match_idx).name);
load(file_path);

end
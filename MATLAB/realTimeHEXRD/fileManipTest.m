filePath = '/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-lineup/';
filePattern = fullfile(filePath, '/*/ff/ff1*.h5'); 
filePattern2 = fullfile(filePath, '/*/ff/ff2*.h5');
files = dir(filePattern);
files2 = dir(filePattern2);
[fileNums] = cellfun(@(x)sscanf(x,'ff1_%d.h5'), {files.name});
[sfNums,sortInd] = sort(fileNums);
sortInd(sfNums<60) = [];
files = files(sortInd);
files2 = files2(sortInd);
{files.name}
{files.folder}

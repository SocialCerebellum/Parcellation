function dwarfs_validate

% create inputbox
answer = inputdlg({'Unique ALE filename string:', 'ALE Folder:', ...
    'ALE p-value (e.g., "pID=..."):', 'ALE minimum Volume:'}, ...
    'Input', 1, ...
    {'WB_Topic #1 = 1 - trait - kopie', 'G:\CerebellumDwarfs\3. ALE_WB', '', ''});
%    {'WB_Topic #', 'G:\CerebellumDwarfs\3. ALE_WB', 'pID=0.0001', '500'});

ALEfilestring = [answer{1} '*.txt'];
ALEfolder = [answer{2} '\'];
ALEpvalue = answer{3};
ALEminVol = answer{4};
ALEpath = [ALEfolder ALEfilestring];

% Select ALE text files

disp('Selecting ALE text files')
% empty ppdirectories
ff = 0;
ffdirectories = {};
% for each filename
TT1 = 0;
% make a list of all ALE files:
ALEfilenames = dir(ALEpath);
ff = (size(ALEfilenames, 1));
for T1 = 1:ff
    try
        % to resolve 'txt' name confusion
        if ~(contains(ALEfilenames(T1).name, 'clust'))
        TT1 = TT1 + 1;
        ffdirectories{TT1}= ALEfilenames(T1).name;
        end
    continue
    end
end
%ff = TT1;
nALEtopics = TT1;
disp(['Number of available ALE topics = ' num2str(nALEtopics)]);

for T2 = 1:nALEtopics
    ALEfile = ffdirectories{T2};
    command = ...
        ['java -cp "C:\Program Files\GingerAle\GingerALE.jar"' ...
        ' org.brainmap.meta.getALE2 "' ALEfolder ALEfile '"' ...
        ' -' ALEpvalue ' -minVol=' ALEminVol ' -mask=MNI_wb_dil.nii'];
    system(command);
    disp(['Completed ALE analysis of file: ' ALEfile])
end
end


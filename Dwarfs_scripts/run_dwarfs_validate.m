function dwarfs_validate
% Reliability analysis
% The data set is randomly split in two parts (termed 'Train' and 'Test', 
% although both are equivalent) with 'Percentage Leave Out' the exact % 
% of studies in the 'Test' part. For each part, an ALE analysis can be 
% conducted. This procedure is repeated for each cross-validation.

% create inputbox
answer = inputdlg({'Unique ALE filename string:', 'Validation Folder', ...
    'Percentage Leave Out (% studies for Test; if 0 "_MNI" file name)' ...
    'Start of Replications (empty for only ALE batch)', 'End of Replications', ...
    'Batch ALE (y or n)'}, ...
    'Input', 1, ...
    {'Topic', 'G:\CerebellumDwarfs\4. CROSSVAL', '50', '14', '20', 'y'});
%    {'Topic #0 ', 'G:\CerebellumDwarfs\4. CROSSVAL', '50', '', '', 'n'});
%    {'Topic #', 'G:\CerebellumDwarfs\4. CROSSVAL', '5', '1', '20', 'y'});

ALEfilestring = [answer{1} '*.txt'];
VALfolder = [answer{2} '\'];
VALleaveP = str2num(answer{3});
VALrepsta = str2num(answer{4});
VALrepend = str2num(answer{5});
ALEbatch = answer{6};

ALEfolder = VALfolder ;
ALEpath = [ALEfolder ALEfilestring];

%======================
% Select ALE text files
%======================

disp('== Select ALE text files ==')
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
        % to resolve name confusions
        if ~(contains(ALEfilenames(T1).name, 'clust')) & ...
           ~(contains(ALEfilenames(T1).name, 'test')) & ...
           ~(contains(ALEfilenames(T1).name, 'ALE')) & ...
           ~(contains(ALEfilenames(T1).name, '_pID')) & ...
           ~(contains(ALEfilenames(T1).name, 'TEMP')) & ...
           ~(contains(ALEfilenames(T1).name, 'Topics on')) & ...
           ~(contains(ALEfilenames(T1).name, 'NoTopic')) & ...
           ~(contains(ALEfilenames(T1).name, '#99')) & ...
           ~(contains(ALEfilenames(T1).name, '_CVAL')) & ...
           ~(contains(ALEfilenames(T1).name, '_MNI')) & ...
           ~(contains(ALEfilenames(T1).name, '%')) & ...
           ~(contains(ALEfilenames(T1).name, 'Clusters'))
        TT1 = TT1 + 1;
        ffdirectories{TT1}= ALEfilenames(T1).name;
        disp(ALEfilenames(T1).name)
        end
    continue
    end
end
ff = TT1;
nALEtopics = TT1;
disp(['Number of available ALE topics = ' num2str(nALEtopics)]);

%===========
% Resampling 
%===========

if VALrepsta > 0
    nsta = VALrepsta;
    nend = VALrepend;
else
    nsta = 1;
    nend= 1;
end

for C1 = nsta:nend
    disp(['============================']);
    disp(['== CROSS-VALIDATION ' num2str(C1) ' ' num2str(VALleaveP) '% ==']);
    disp(['============================']);

    %=======================================================
    % Read ALE txt files and leave-%-studies-out for testing
    %=======================================================

    for T2 = 1:nALEtopics
        % read files
        ALEtextfile = ffdirectories{T2};
        copyfile([VALfolder ALEtextfile], [VALfolder 'TEMP.txt']);
        file = fopen([VALfolder 'TEMP.txt']);
        trai = fopen([VALfolder 'TEMPtrain.txt'], 'w');
        test = fopen([VALfolder 'TEMPtest.txt'], 'w');

        % defaults
        nStudiesTest = 0;
        nStudiesTrain = 0;
        nStudies = 0;

        Line = fgets(file);       
        while ischar(Line)
            Line = fgets(file);
            if length(Line) > 1
                if contains(Line, '//') & contains(Line, 'Subjects=')
                nStudies = nStudies + 1;
                end
            end
        end

        % read first line
        frewind(file);
        Line = fgets(file);
        
        % add reference system if absent
        if ~contains(Line, 'Reference=MNI')
            Line0 = '// Reference=MNI';
            fprintf(trai, '%s\n', Line0);
            fprintf(test, '%s\n', Line0);
        end
        % remove 'Topic' line
        if contains(Line, 'Topic #')
            Line0 = ' ';
            fprintf(trai, '%s\n', Line0);
            fprintf(test, '%s\n', Line0);
        end

        while ischar(Line)
            Line_prev = Line;
            Line = fgets(file);
            if length(Line) > 1
                % seek study line ('//')
                if contains(Line, '//') & contains(Line, 'Subjects=')
                    % random test or train
                    if unidrnd(100) < VALleaveP
                        ff = test;
                        nStudiesTest = nStudiesTest + 1;
                    else
                        ff = trai;
                        nStudiesTrain = nStudiesTrain + 1;
                    end
                    % set limits on test or train
                    if nStudiesTest  > (VALleaveP / 100 * nStudies) + .5 
                        ff = trai;
                        nStudiesTest = nStudiesTest - 1;
                        nStudiesTrain = nStudiesTrain + 1;
                    end
                    if nStudiesTrain > ((100 - VALleaveP) / 100 * nStudies) + .5
                        ff = test;
                        nStudiesTrain = nStudiesTrain - 1;
                        nStudiesTest = nStudiesTest + 1;
                    end
                    % blank, author & subject lines  
                    fprintf(ff, '%s\n', ' '); 
                    fprintf(ff, '%s', Line_prev); 
                    fprintf(ff, '%s', Line);
                else
                    % coordinate lines 
                    if ~contains(Line, 'Topic #') & ~contains(Line, '//')
                        if ff == test | ff == trai
                            fprintf(ff, '%s', Line);
                        end
                    end
                end
            end
        end 
        % last line
        if length(Line) > 1
            fprintf(ff, '%s', Line);
        end
        fclose(file);
        fclose(trai);
        fclose(test);

        ALEtextpart = ALEtextfile(1:length(ALEtextfile) - 4);
        disp([ALEtextpart ' - Studies for training = ' num2str(nStudiesTrain) ...
              ' for testing = ' num2str(nStudiesTest)])
        if VALleaveP > 0
            copyfile([VALfolder 'TEMPtrain.txt'], [VALfolder ALEtextpart '_Ctrain' num2str(C1) '_' num2str(VALleaveP) '%.txt']);
            copyfile([VALfolder 'TEMPtest.txt'], [VALfolder ALEtextpart '_Ctest' num2str(C1) '_' num2str(VALleaveP) '%.txt']);
        else
            copyfile([VALfolder 'TEMPtrain.txt'], [VALfolder ALEtextpart '_MNI.txt']);
        end

        %======================================
        % Leave-%-out cross-validation by ALE
        %======================================

        if ~contains(ALEbatch, 'n') 
            % ALE of train & test
            for i = 1:2
                ALEtextpart = ALEtextfile;
                ALEtextpart = ALEtextpart(1:length(ALEtextpart) - 4); 
                if VALleaveP > 0
                   if i == 1; ALEtextpart = [ALEtextpart  '_Ctrain' num2str(C1)]; end
                   if i == 2; ALEtextpart = [ALEtextpart  '_Ctest' num2str(C1)]; end
                   ALEfilepart = [ALEtextpart '_' num2str(VALleaveP) '%.txt'];
                else
                   ALEfilepart = [ALEtextpart '_MNI.txt'];
                end

                % ALE of train & test | ALE only once for MNI.txt
                if VALleaveP > 0 | VALleaveP <= 0 & i == 1
                    disp([' ... Starting ALE cross-validation ' num2str(C1) ' of file: ' ALEtextpart])
                    command = ...
                        ['java -cp "C:\Program Files\GingerAle\GingerALE.jar"' ...
                        ' org.brainmap.meta.getALE2 "' VALfolder ALEfilepart '"' ...
                        ' -mask=MNI_wb_dil.nii'];
        %                ' -mask=MNI_wb_dil.nii -pID=0.0001 -minVol=500'];
                    system(command);
                    disp([' ... Completed ALE cross-validation ' num2str(C1) ' of file: ' ALEtextpart])
                end
            end

        end
    end
end
fclose('all');
delete([VALfolder 'TEMPtrain.txt']);
delete([VALfolder 'TEMPtest.txt']);
delete([VALfolder 'TEMP.txt']);
end
    

function ALEInput

% create inputbox
answer = inputdlg({'Excel Database:', 'Excel Database Folder:', ...
                   'ALE File Prefix (CE_, WB_):', 'ALE Folder:', ...
                   'Topics (y, o[pposite], c[onjoint], a[ll: y + o + c])'}, ...
    'Input', 1, ...
    {'database_CE_Valid.xlsx', 'G:\CerebellumDwarfs\1. DATABASE', ...
    'CE_', 'G:\CerebellumDwarfs\3. ALE_2', 'a'});
%    {'database_WB_Valid.xlsx', 'G:\CerebellumDwarfs\1. DATABASE', ...
%    'WB_', 'G:\CerebellumDwarfs\3. ALE_2', 'a'});

wb_read = answer{1};
folder = [answer{2} '\'];
prefix = answer{3};
ALEfolder = [answer{4} '\'];
Topics = answer{5};
nTopics = 22;

mkdir(ALEfolder);

disp(['reading database: ' folder wb_read])
[num, txt, raw] = xlsread([folder wb_read]);
%[num, txt, raw] = xlsread([folder wb_read], 'A1:BZ500');
disp(['reading database: ' folder wb_read ' ended'])

%%Row2Write = 25;
Row2Start = 5;
Row2End = size(num, 1);
Column2Read = 26 + 1 - 6;  

% ===================================================
% ALE input for each topic - Loop through all nTopics 
% ===================================================
if Topics == 'y' | Topics == 'a' 

    for NextTopic = 1:nTopics
        Column2Read = Column2Read + 6;
        prior_article = 0;
        nArticles = 0;
        nSubjects = 0;
        nFoci = 0;

        % For each topic
        sTopic = num2str(num(1, Column2Read));
        sTerm1 = char(txt(1, Column2Read + 1));
        sNextTopic = num2str(NextTopic);

        ff = [ALEfolder prefix 'Topic #' sNextTopic ' = ' sTopic ' - ' sTerm1 '.txt'];
        fileID = fopen(ff, 'w');

        fprintf(fileID,'// Reference=MNI');
        fprintf(fileID,'\n');
        fprintf(fileID,['// Topic #' sNextTopic ' = ' sTopic ' - ' sTerm1]);
        fprintf(fileID,'\n');

        % Check if coordinates are valid for given terms
        for NextRow = Row2Start:Row2End
            if num(NextRow, Column2Read) == 1 

                % Identify a new study
                if num(NextRow, 1) ~= prior_article 
                    sAuthors = char(txt(NextRow, 11));
                    sYear = num2str(num(NextRow, 12));
                    sSubjects = num2str(num(NextRow, 15));

                    fprintf(fileID,'\n');
                    fprintf(fileID,['// ' sAuthors ' (' sYear ')' '\n']);
                    fprintf(fileID,['// Subjects=' sSubjects '\n']);

                    prior_article = num(NextRow, 1);              
                    nArticles = nArticles + 1;
                    nSubjects = nSubjects + str2num(sSubjects);
                end 
                Coord_x = num(NextRow, 3);
                Coord_y = num(NextRow, 4);
                Coord_z = num(NextRow, 5);
                C = [Coord_x, Coord_y, Coord_z];
                fprintf(fileID, '%4i', C);
                fprintf(fileID,'\n');
                nFoci = nFoci + 1;
            end 
        end
        disp([ff ' with Articles = ' num2str(nArticles) ' Foci = ' num2str(nFoci) ' Subjects = ' num2str(nSubjects)])
        fclose(fileID);
    end
end
    
% ==================================================
% ALE input without topic - Loop through all nTopics 
% ==================================================
if Topics == 'o' | Topics == 'a'

    Column2Read = 26 + 1 - 6;  
    for NextTopic = 1:nTopics
        Column2Read = Column2Read + 6;
        prior_article = 0;
        nArticles = 0;
        nSubjects = 0;
        nFoci = 0;

        % For each topic
        sTopic = num2str(num(1, Column2Read));
        sTerm1 = char(txt(1, Column2Read + 1));
        sNextTopic = num2str(NextTopic);

        ff = [ALEfolder prefix 'NoTopic #' sNextTopic ' = ' sTopic ' - ' sTerm1 '.txt'];
        fileID = fopen(ff, 'w');

        fprintf(fileID,'// Reference=MNI');
        fprintf(fileID,'\n');
        fprintf(fileID,['// NoTopic #' sNextTopic ' = ' sTopic ' - ' sTerm1]);
        fprintf(fileID,'\n');

        % For all rows
        for NextRow = Row2Start:Row2End
            ifValid = 0;
            ValidColumn2Read = 26 + 1 - 6;  
            % Check if coordinates are valid for at least one term
            for NextValidTopic = 1:nTopics
                ValidColumn2Read = ValidColumn2Read + 6;
                if num(NextRow, ValidColumn2Read) == 1 
                    ifValid = 1;
                end
            end
            % Coordinates are valid for at least one term & NOT for given term
            if ifValid == 1 & num(NextRow, Column2Read) ~= 1

                % Identify a new study
                if num(NextRow, 1) ~= prior_article 
                    sAuthors = char(txt(NextRow, 11));
                    sYear = num2str(num(NextRow, 12));
                    sSubjects = num2str(num(NextRow, 15));

                    fprintf(fileID,'\n');
                    fprintf(fileID,['// ' sAuthors ' (' sYear ')' '\n']);
                    fprintf(fileID,['// Subjects=' sSubjects '\n']);

                    prior_article = num(NextRow, 1);              
                    nArticles = nArticles + 1;
                    nSubjects = nSubjects + str2num(sSubjects);
                end 
                Coord_x = num(NextRow, 3);
                Coord_y = num(NextRow, 4);
                Coord_z = num(NextRow, 5);
                C = [Coord_x, Coord_y, Coord_z];
                fprintf(fileID, '%4i', C);
                fprintf(fileID,'\n');
                nFoci = nFoci + 1;
            end
        end
        disp([ff ' with Articles = ' num2str(nArticles) ' Foci = ' num2str(nFoci) ' Subjects = ' num2str(nSubjects)])
        fclose(fileID);
    end
end
    
% ====================================================
% ALE input with all topics - Loop through all nTopics 
% ====================================================
if Topics == 'c' | Topics == 'a'

    ff = [ALEfolder prefix 'Topic #0 = all.txt'];
    fileID = fopen(ff, 'w');

    fprintf(fileID,'// Reference=MNI');
    fprintf(fileID,'\n');
    fprintf(fileID,'// All topics');
    fprintf(fileID,'\n');
    
    prior_article = 0;
    nArticles = 0;
    nSubjects = 0;
    nFoci = 0;
    
    % For all rows
    for NextRow = Row2Start:Row2End
        ifValid = 0;
        ValidColumn2Read = 26 + 1 - 6;  
        % Check if coordinates are valid for at least one term
        for NextValidTopic = 1:nTopics
            ValidColumn2Read = ValidColumn2Read + 6;
            if num(NextRow, ValidColumn2Read) == 1 
                ifValid = 1;
            end
        end
        % Coordinates are valid for at least one term
        if ifValid == 1
            
            % Identify a new study
            if num(NextRow, 1) ~= prior_article 
                sAuthors = char(txt(NextRow, 11));
                sYear = num2str(num(NextRow, 12));
                sSubjects = num2str(num(NextRow, 15));

                fprintf(fileID,'\n');
                fprintf(fileID,['// ' sAuthors ' (' sYear ')' '\n']);
                fprintf(fileID,['// Subjects=' sSubjects '\n']);

                prior_article = num(NextRow, 1);              
                nArticles = nArticles + 1;
                nSubjects = nSubjects + str2num(sSubjects);
            end 
            Coord_x = num(NextRow, 3);
            Coord_y = num(NextRow, 4);
            Coord_z = num(NextRow, 5);
            C = [Coord_x, Coord_y, Coord_z];
            fprintf(fileID, '%4i', C);
            fprintf(fileID,'\n');
            nFoci = nFoci + 1;
        end
    end
    disp([ff ' with Articles = ' num2str(nArticles) ' Foci = ' num2str(nFoci) ' Subjects = ' num2str(nSubjects)])
    fclose(fileID);
end 

end


function dwarfs_DCBC
% run DCBC analysis

%
% ENTER TITLE, CLUSTERS AND THEIR LABELS HERE BELOW
%

Title =    'Ward 61 on MDTB (euclidean)';													
Clusters = [20	54	0	0	0	0	0	;	...
			22	56	29	61	0	0	0	;	...
			43	44	55	0	0	0	0	;	...
			4	35	15	25	57	33	48	;	...
			3	34	19	0	0	0	0	;	...
			26	58	27	28	59	60	0	;	...
			45	47	46	51	0	0	0	;	...
			1	30	9	17	50	0	0	;	...
			2	36	13	7	31	40	0	;	...
			5	10	18	0	0	0	0	;	...
			6	16	14	21	0	0	0	;	...
			8	38	39	0	0	0	0	;	...
			32	49	37	0	0	0	0	;	...
			11	42	12	24	41	0	0	;	...
			23	52	53	0	0	0	0	];	
Labels   = {'empathy              ',	...
			'word                 ',	...
			'prediction           ',	...
			'rest-video           ',	...
			'finger-sequence      ',	...
			'visual-search        ',	...
			'response-alternatives',	...
			'action-observation   ',	...
			'movies               ',	...
			'conflict             ',	...
			'emotion              ',	...
			'math-rotation        ',	...
			'biological-motion    ',	...
			'memory               ',	...
			'spatial              '};	

%Title =    'CE+WB Complete (#98) on Z (euclidean)';							
%Clusters = [1	5  15	0 0;   2  3  18 19 0;  6   7  16  20 13; 11  0	0	0 0;   ...
%            4   0	0	0 0;  10  0	 0	 0 0;  8	0	0	0 0; 14	21	0	0 0;  9	0	0	0 0; ...
%            22 12  17   0 0];	
%Labels   = {'empathy',      'conflict',     'decision',       'face',    ...
%            'emotion',      'objects',      'encoding',       'attention',     'word', ...
%            'sensorimotor'};	
%
%Title =    'CE+WB Complete (#98) on Z (euclidean) 7-networks';							
%Clusters = [1	5  15	0 0 0 0;   2  3 18 19 0 0 0;  6 7 16 20 13 10 8;   ...
%            4  11	0	0 0 0 0;  14 21	 0  0 0 0 0;  9	0  0  0  0  0 0; ...
%            22 12  17   0 0 0 0];	
%Labels   = {'empathy',         'conflict',         'decisionXL', ...
%            'emotionXL',       'attention',        'word', ...
%            'sensorimotor'};	

% create inputbox
answer = inputdlg({'Unique Topics filename string (before .nii):', ...
    'Topics folder:', ...
    'Create Volume Flatmaps (y, n, o[nly]):', ...
    'Colormap Filename (in ALE folder):', ...
    'Cluster Title:', ...
    'EquiDistance start value:', ...
    'EquiDistance end value:'}, ...
    'Input', 1, ...
    {'_Z', 'G:\CerebellumDwarfs\3. ALE', 'n', 'Buckner-likeColors.txt', Title, '5', '80'});
%    {'Set', 'G:\CerebellumDwarfs\MDTB', 'n', 'Buckner-likeColors-MDTB.txt', Title, '5', '80'});
%    {'_CVAL1_5%_Z', 'G:\CerebellumDwarfs\4. CROSSVAL', 'n', 'Buckner-likeColors.txt', Title});
%    {'_Z', 'G:\CerebellumDwarfs\3. ALE_WB', 'n', 'Buckner-likeColors.txt', Title});
%    {'Set', 'G:\CerebellumDwarfs\MDTB', 'n', 'MDTB-10Regions.txt', Title});

if size(answer{1}, 1) ~= 0
    ALEfilestring = ['*' answer{1} '*.nii'];
else
	ALEfilestring = ['*.nii']; 
end
ALEfolder = [answer{2} '\'];
Flatmaps = answer{3};
Colormap = answer{4};
Title = answer{5};
EquiDistanceStart = str2num(answer{6});
EquiDistanceEnd = str2num(answer{7});

% defaults
nALEtopics = 0;
ALEpath = [ALEfolder ALEfilestring];
Colormap = [ALEfolder Colormap];
nClusters = size(Clusters, 1);

%=========
% Colormap
%=========

[CM(:, 1) CM(:, 2) CM(:, 3) CM(:, 4) CMlabels] = readvars(Colormap);
Colors = [CM(:,2) CM(:,3) CM(:,4)]/255;
%Topics = char(CMlabels);
Topics = convertCharsToStrings(CMlabels);

% Reordering of cluster colors for 'Buckner-likeColors.txt'
if contains(answer{4}, 'Buckner-likeColors')
    CC = [];
    for i = 1:size(Topics, 1)
        for j = 1:size(Labels, 2)
            s1 = Topics(i, :);
            s2 = strip(Labels(j));
            if contains(s1, s2)
                CC(j,:) = Colors(i,:); 
            end
        end
    end
    Colors = CC;
end
Topics = Labels';

for i = 1:size(Colors, 1)
    Colors2Map(i, :) = i;
    Colors2Map(i, [2:4]) = Colors(i, [1:3]) * 255;
end
save([ALEfolder 'TEMPcolormap.txt'], 'Colors2Map', '-ascii');
Colormap = [ALEfolder 'TEMPcolormap.txt'];

% activate SUIT toolbox
if Flatmaps ~= 'n'
    spm_suit
end

%=============
% Select Files
%=============

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
        if ~(contains(ALEfilenames(T1).name, 'ALE')) & ...
           ~(contains(ALEfilenames(T1).name, '_pID')) & ...
           ~(contains(ALEfilenames(T1).name, 'TEMP')) & ...
           ~(contains(ALEfilenames(T1).name, 'Topics on')) & ...
           ~(contains(ALEfilenames(T1).name, 'NoTopic')) & ...
           ~(contains(ALEfilenames(T1).name, '#0')) & ...
           ~(contains(ALEfilenames(T1).name, '#99')) & ...
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
disp(['Number of available topics = ' num2str(nALEtopics)]);

%===============================
% Assign clusters & Read Volumes
%===============================

% cluster assignments
for V1 = 1:size(Clusters, 1)
    for j = 1:size(Clusters, 2)
        if (Clusters(V1, j)) > 0
        ClusterAssignment(Clusters(V1, j)) = V1;
        end
    end
end
    
% read volumes
X = [];
for T2 = 1:nALEtopics
    ALEfilename = ffdirectories{T2};
    info = niftiinfo([ALEfolder ALEfilename]);
    V = niftiread(info);
    % write in linear (volume) matrix X
    V1 = 0;
    for x = 1:size(V, 1)
        for y = 1:size(V, 2)
            for z = 1:size(V, 3)
                V1 = V1 + 1;
                X(T2, V1) = V(x, y, z);
            end
        end
    end
    
%================
% Create Flatmaps
%================
    
    if Flatmaps ~= 'n'
        Data = suit_map2surf([ALEfolder ALEfilename], 'space','SPM','stats',@minORmax);
        disp(['Create Flatmap ' ALEfilename ' = ' num2str(T2)]) 
        figure();
        suit_plotflatmap(Data, 'cscale', [], 'cmap', hot);
        ALEfile = ALEfilename(1:size(ALEfilename, 2) - 4);
        title(ALEfile);
        exportgraphics(gcf, [ALEfolder 'TT Flatmap' ALEfile '.png'])
    end
end
if Flatmaps == 'o'
    return
end

%==========================
% Rescale & Winner-take-all
%==========================

% rescale volumes
for T3 = 1:nALEtopics
    Xmax(T3) = 0;
    Xmin(T3) = 0;
    for V1 = 1:size(X, 2)
        if Xmax(T3) < X(T3, V1)
            Xmax(T3) = X(T3, V1);
        end
        if Xmin(T3) > X(T3, V1)
            Xmin(T3) = X(T3, V1);
        end
    end
    for V1 = 1:size(X, 2)
        % normalize
        X(T3, V1) = X(T3, V1) / Xmax(T3);
        % standardize
        %X(T3, V1) = X(T3, V1) + Xmin(T3);
        %X(T3, V1) = X(T3, V1) / (Xmax(T3) + Xmin(T3));
    end
end
Title = [Title ' Rescaled '];
% winner-take-all
for V1 = 1:size(X, 2)
    XmaxTopic(V1) = 0;
    XmaxCluster(V1) = 0; 
    maxy = 0;
    for T2 = 1:size(X, 1)
        if maxy < X(T2, V1)
            maxy = X(T2, V1);
            XmaxTopic(V1) = T2;
            XmaxCluster(V1) = ClusterAssignment(T2);
        end
    end
end
disp(['Number of clusters = ' num2str(nClusters)])

%======================================
% Task activation used to see if there
% is activation for a given voxel)
%======================================
   
Xactivation = [];
for V1 = 1:size(X, 2)
    Xactivation(V1) = 0;
    for T4 = 1:nALEtopics
        Xactivation(V1) = Xactivation(V1) + X(T4, V1);
    end
end

%=================================
% Compare voxels at given distance
%=================================

% defautls for Excel
xlsfilename = [ALEfolder 'DCBC_' Title '.xls'];
sheet = 1;
r = 1;
        
% define equidistance for bins
for EquiDistance = EquiDistanceStart:5:EquiDistanceEnd
    r = r + 2;
    % get dimensionality 
    ALEfilename = ffdirectories{1};
    info = niftiinfo([ALEfolder ALEfilename]);
    V = niftiread(info);
    % reading voxel size 
    Vmm = info.PixelDimensions(1);
    EquiDistance = round(EquiDistance/Vmm);
    disp(['EquiDistance = ' num2str(EquiDistance * Vmm)])    
    
    % defaults
    XtargetW = [];
    XtargetB = [];
    XcomparW = [];
    XcomparB = [];
    XtargetWC = [];
    XcomparWC = [];
    j = 0;
    k = 0;
    V1 = 0;
    progress = 0;
    progress_old = 0;
    
    % target voxels in V 
    % write in linear (volume) matrix X
    disp('Determining equidistant bins ... ')    
    for x1 = 1:size(V, 1)
        for y1 = 1:size(V, 2)
            for z1 = 1:size(V, 3)
                V1 = V1 + 1;
                if Xactivation(V1) > 0

                    % comparison voxels in V 
                    % write in linear (volume) matrix X
                    V2 = 0;
                    for x2 = 1:size(V, 1)
                        for y2 = 1:size(V, 2)
                            for z2 = 1:size(V, 3)
                                V2 = V2 + 1;
                                % stay in equidistant box (this is the
                                % fastest if condition)
                                if x2 > x1 - EquiDistance - Vmm & x2 < x1 + EquiDistance + Vmm & ...
                                   y2 > y1 - EquiDistance - Vmm & y2 < y1 + EquiDistance + Vmm & ...
                                   z2 > z1 - EquiDistance - Vmm & z2 < z1 + EquiDistance + Vmm & ...
                                   Xactivation(V2) > 0

                                    % compute distance    
                                    euclidian = realsqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
                                    %if euclidian > EquiDistance - 0.02 & euclidian < EquiDistance + 0.02
                                    if euclidian == EquiDistance
                                        %%%disp(['linear voxels = ' num2str(V1) ' & ' num2str(V2) ' with parcellations ' num2str(XmaxCluster(V1)) ' & ' num2str(XmaxCluster(V2))])
                                        % stack activations in equidistant bins for all topics 
                                        for T4 = 1:nALEtopics
                                            if XmaxCluster(V1) == XmaxCluster(V2)
                                                % within parcellations  
                                                j = j + 1;
                                                % stacked per parcellation
                                                XtargetWC(XmaxCluster(V1),j) = X(T4, V1);
                                                XcomparWC(XmaxCluster(V2),j) = X(T4, V2);
                                                %%%disp(['topic = ' num2str(T4) ' with W activation ' num2str(XtargetWC(XmaxCluster(V1),j)) ' & ' num2str(XcomparWC(XmaxCluster(V2),j))  ])
                                                % stacked overall
                                                XtargetW(j) = X(T4, V1);
                                                XcomparW(j) = X(T4, V2);
                                                %disp(['topic = ' num2str(T4) ' with W activation ' num2str(XtargetW(j)) ' & ' num2str(XcomparW(j))  ])
                                            else
                                                % between parcellations                                                   
                                                k = k + 1;
                                                % stacked overall
                                                XtargetB(k) = X(T4, V1);
                                                XcomparB(k) = X(T4, V2);
                                                %disp(['topic = ' num2str(T4) ' with B activation ' num2str(XtargetB(k)) ' & ' num2str(XcomparB(k))  ])
                                            end
                                        end
                                        
                                        progress = round(x1 / size(V, 1) * 100);
                                        if progress_old < progress
                                        %if mod(j, size(V, 1) * size(V, 2) * size(V, 3) / 10) == 0 
                                            disp(['   ... ' num2str(progress) ' %'])
                                            progress_old = progress;
                                        end 
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % create pairwise missing values
    for T4 = 1:size(XtargetWC, 1)
        for kk = 1:size(XtargetWC, 2)
            if XtargetWC(T4, kk) == 0; XtargetWC(T4, kk) = NaN; end
            if XcomparWC(T4, kk) == 0; XcomparWC(T4, kk) = NaN; end
        end
    end
                                        
    % Write on Excel
    disp (['Equidistance = ' num2str(EquiDistance * Vmm)])

    [rWithin,  pWithin]  = corrcoef(XtargetW, XcomparW);
    [rBetween, pBetween] = corrcoef(XtargetB, XcomparB);
    
    xlswrite(xlsfilename, {['Equidistance = ' num2str(EquiDistance * Vmm)]}, sheet, [char(64 + 1)  num2str(r)]);
    r = r + 1;
    xlswrite(xlsfilename, {'Overall Within'}, sheet, [char(64 + 1)  num2str(r)]);
    xlswrite(xlsfilename, rWithin(2, 1), sheet, [char(64 + 3)  num2str(r)]);
    xlswrite(xlsfilename, pWithin(2, 1), sheet, [char(64 + 4)  num2str(r)]);
    xlswrite(xlsfilename, {'Overall Between'}, sheet, [char(64 + 6)  num2str(r)]);
    xlswrite(xlsfilename, rBetween(2, 1), sheet, [char(64 + 8)  num2str(r)]);
    xlswrite(xlsfilename, pBetween(2, 1), sheet, [char(64 + 9)  num2str(r)]);
       
    %for i = 1:nClusters
    for i = 1:size(XtargetWC, 1)
        disp(['Cluster/Parcellation = ' num2str(i)])
        [rWithin,  pWithin]  = corrcoef(XtargetWC(i,:), XcomparWC(i,:), 'Rows', 'pairwise');
        r = r + 1;
        xlswrite(xlsfilename, {'Parcellation Within '}, sheet, [char(64 + 1)  num2str(r)]);
        xlswrite(xlsfilename, Labels(i), sheet, [char(64 + 2)  num2str(r)]);
        xlswrite(xlsfilename, rWithin(2, 1), sheet, [char(64 + 3)  num2str(r)]);
        xlswrite(xlsfilename, pWithin(2, 1), sheet, [char(64 + 4)  num2str(r)]);
    end

    clear XtargetW;
    clear XtargetB;
    clear XcomparW;
    clear XcomparB;
end 

end
   
    

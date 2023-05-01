function dwarfs_winner
%
% ENTER TITLE, CLUSTERS AND THEIR LABELS HERE BELOW
%

Title =    'Complete+Ward on Z (euclidean)';
Clusters = [ 1  5 15  0; 18 19  0  0;  2  3  0  0;  6  7  0  0; 12 17  0  0; 16 20 0 0; ...
             4 11  0  0; 10 13  0  0;  8  0  0  0; 14 21  0  0;  9  0  0  0 ];           
Labels   = {'empathy', 'fear',       'conflict',  'decision',   'training', 'working',...
            'emotion', 'objects',    'encoding',  'attention',  'word',   };

Title =    'WB Ward+Complete on Z (euclidean)';
Clusters = [ 4	11	0	0;	1	18	13	0;	5	15	0	0;	2	7	20	0;	 6   0   0   0;...
            19	0	0	0;	3	0	0	0;	8	10	0	0;	14	16	0	0; ...
             9	0	0	0;	12	0	0	0;	17	0	0	0;	21	0	0	0];
Labels   = {'emotion'	 ,'fear'		 ,'empathy'      ,'conflict'     ,'decision' ...
            'inhibition' ,'location'     ,'encoding'     ,'attention'    , ...
            'word'       ,'reward'       ,'training'     ,'action'       };			

Title =    'Complete+Ward on Z (euclidean)';
Clusters = [ 1  5 15  0; 18 19  0  0;  2  3  0  0;  6  7  0  0; 12 17  0  0; 16 20 0 0; ...
             4 11  0  0; 10 13  0  0;  8  0  0  0; 14 21  0  0;  9  0  0  0 ];           
Labels   = {'empathy', 'fear',       'conflict',  'decision',   'training', 'working',...
            'emotion', 'objects',    'encoding',  'attention',  'word'   };

Title =    '7-network of Complete on Z (euclidean)';
Clusters = [ 1  5 15  0  0  0  0  0  0; 18 19  4 11 0 0 0 0 0;  2  3  0  0 0 0 0 0 0;...
             6  7 12 17 16 20 10 13  8; 14 21  0  0 0 0 0 0 0;  9  0  0  0 0 0 0 0 0];           
Labels   = {'empathy',                  'fear',                'conflict', ...
            'decision',                 'attention',           'word'   };

Title =    'Complete on Z (euclidean)';
Clusters = [ 1  5 15  0; 18 19  0  0;  2  3  0  0;  6  7  0  0; 12 17 16 20; ...
             4 11  0  0; 10 13  0  0;  8  0  0  0; 14 21  0  0;  9  0  0  0 ];           
Labels   = {'empathy', 'fear',       'conflict',  'decision',   'working',...
            'emotion', 'objects',    'encoding',  'attention',  'word'   };

Title =    'Ward 64 on MDTB (euclidean)';														
Clusters = [25	59	32	64	23	57	0	0	0	0	;	...
			4	33	12	20	53	0	0	0	0	0	;	...
			5	39	16	34	43	9	19	17	24	21	;	...
			3	29	61	30	31	62	63	0	0	0	;	...
			7	18	38	28	60	36	51	0	0	0	;	...
			2	10	46	47	58	8	13	0	0	0	;	...
			11	41	42	35	52	40	0	0	0	0	;	...
			14	45	26	15	27	44	0	0	0	0	;	...
			1	6	37	22	0	0	0	0	0	0	;	...
			48	50	49	54	55	56	0	0	0	0	];  	 									
Labels   = {'empathy',					...
        	'motion-perception',		...
        	'scenes',					...
        	'visual',					...
        	'motor',					...
        	'prediction',				...
        	'spatial-mental',			...
        	'working',					...
        	'motor',					...
        	'spatial-search'};        

Title =    'WB Ward+Complete on Z (euclidean)';
Clusters = [ 4	11	0	0;	1	18	13	0;	5	15	0	0;	2	7	20	0;	 6   0   0   0;...
            19	0	0	0;	3	0	0	0;	8	10	0	0;	14	16	0	0; ...
             9	0	0	0;	12	0	0	0;	17	0	0	0;	21	0	0	0];
Labels   = {'emotion'	 ,'fear'		 ,'empathy'      ,'conflict'     ,'decision' ...
            'inhibition' ,'location'     ,'encoding'     ,'attention'    , ...
            'word'       ,'reward'       ,'training'     ,'action'       };			

%Title =    'CE+WB 11 Ward on Z (euclidean)';							
%Clusters = [1	5	15	0;	2	3	0	0; 18	19	0	0; 	4	0	0	0; ...
%    		11	0	0	0;  10  13	0	0;  8	0	0	0;  9	0	0	0; ...
%    		14	21	0	0;  6	7	16 20; 12	17	0	0];	
%Labels   = {'empathy'   ,	'conflict',     'fear', 	'emotion',	...
%            'face'	,       'objects',  	'encoding',	'word'	,	...
%        	'attention'	,	'decision', 	'reward'	};	

Title =    'CE+WB Ward on Z (euclidean)';							
Clusters = [1	5	15	0;	18	19	0	0; 	2	3	0	0;  6	7	16 20; 12	17	0	0; ...
            4  11	0	0;  10  13	0	0;  8	0	0	0; 14	21	0	0;  9	0	0	0; ...
            22  0   0   0];	
Labels   = {'empathy',      'fear', 	'conflict',     'decision',      'reward', ...
            'emotion',      'objects',  'encoding',     'attention',     'word', ...
            'sensorimotor'};	

Title =    'CE+WB Ward (10) on Z (euclidean)';							
Clusters = [1	5	15	0;	2	3  18  19;  6	7	16 20; 12	17	0	0; ...
            4  11	0	0;  10  13	0	0;  8	0	0	0; 14	21	0	0;  9	0	0	0; ...
            22  0   0   0];	
Labels   = {'empathy',      'conflict',   'decision',     'reward', ...
            'emotion',      'objects',    'encoding',     'attention',     'word', ...
            'sensorimotor'};	

Title =    'CE+WB Ward (#98) on Z (euclidean)';							
Clusters = [1	5  15	0;	2	3  18  19;  6	7  16  20;  ...
            4  11	0	0;  10  0	0	0;  13  0	0	0;  8	0	0	0; 14	21	0	0;  9	0	0	0; ...
            22 12  17  0];	
Labels   = {'empathy',      'conflict',   'decision',      ...
            'emotion',      'objects',    'events',      'encoding',     'attention',     'word', ...
            'sensorimotor'};	

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

% create inputbox
answer = inputdlg({'Unique Topics filename string (before .nii):', ...
    'Topics folder:', ...
    'Create Volume Flatmaps (y, n, o[nly]):', ...
    'Colormap Filename (in ALE folder):', ...
    'Figure (t = topics, c = clusters):', ...
    'Smoothing filter (1, 2, ... or empty):', ...
    'Cluster Title:'}, ...
    'Input', 1, ...
    {'Set', 'G:\CerebellumDwarfs\MDTB', 'n', 'Buckner-likeColors-MDTB.txt', 'c', '2', Title});
%    {'_CVAL1_5%_Z', 'G:\CerebellumDwarfs\4. CROSSVAL', 'n', 'Buckner-likeColors.txt', 'c', '', Title});
%    {'_Z', 'G:\CerebellumDwarfs\3. ALE', 'n', 'Buckner-likeColors.txt', 'c', '', Title});
%    {'_Z', 'G:\CerebellumDwarfs\3. ALE_WB', 'n', 'Buckner-likeColors.txt', 'c', '', Title});
%    {'Set', 'G:\CerebellumDwarfs\MDTB', 'n', 'MDTB-10Regions.txt', 'c', '', Title});

if size(answer{1}, 1) ~= 0
    ALEfilestring = ['*' answer{1} '*.nii'];
else
	ALEfilestring = ['*.nii']; 
end
ALEfolder = [answer{2} '\'];
Flatmaps = answer{3};
Colormap = answer{4};
Graph = answer{5};
Smoothing = answer{6};
Title = answer{7};

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
if Graph == 'c' & contains(answer{4}, 'Buckner-likeColors')
    CC = [];
    for i = 1:size(Topics, 1)
        for j = 1:size(Labels, 2)
            s1 = Topics(i, :);
            %s2 = char(Labels(j));
            s2 = strip(Labels(j));
            if contains(s1, s2)
                CC(j,:) = Colors(i,:); 
            end
        end
    end
    Colors = CC;
end
if Graph == 'c'
    Topics = Labels';
end;

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
    
% activate SUIT toolbox
if Flatmaps ~= 'n'
    spm_suit
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
        disp(['Create Flatmap ' ALEfilename ' = ' num2str(T2)]) 
        figure();
        Data = suit_map2surf([ALEfolder ALEfilename], 'space','SPM','stats',@minORmax);
        suit_plotflatmap(Data, 'cscale', [], 'cmap', hot);
        ALEfile = ALEfilename(1:size(ALEfilename, 2) - 4);
        title(ALEfile);
        exportgraphics(gcf, [ALEfolder 'Flatmap' ALEfile '.png'])
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
    max = 0;
    for T2 = 1:size(X, 1)
            if max < X(T2, V1)
                max = X(T2, V1);
                XmaxTopic(V1) = T2;
                XmaxCluster(V1) = ClusterAssignment(T2);
            end
    end
end
disp(['Number of clusters = ' num2str(nClusters)])

%================
% Output Graphics
%================

if Flatmaps == 'n'
    %spm_suit
end
ALEfile = [' on ' answer{1}(1:size(answer{1}, 2))];

if Graph == 't'; g = 1; end
if Graph == 'c'; g = 2; end

for V1 = g:g
    % name of graph
    if V1 == 1; GraphName = 'Topics'; end
    if V1 == 2; GraphName = 'Clusters'; end
    % rewrite from linear X to volume matrix V
    V = niftiread(info);
    i = 0;
    for x = 1:size(V, 1)
        for y = 1:size(V, 2)
            for z = 1:size(V, 3)
                i = i + 1;
                if V1 == 1; V(x, y, z) = single(XmaxTopic(i)); end
                if V1 == 2; V(x, y, z) = single(XmaxCluster(i)); end
            end
        end
    end
    
    %Vtest= V(20:40, 20:35, 40)    
    if ~isempty(Smoothing)
    s = str2num(Smoothing);
    VV = V;
    % smoothing SS times
    for SS = 1:10
        rx = randperm(size(V, 1) - 2*s);
        ry = randperm(size(V, 2) - 2*s);
        rz = randperm(size(V, 3) - 2*s);
        for x = 1:size(V, 1) - 2*s
            for y = 1:size(V, 2) - 2*s
                for z = 1:size(V, 3) - 2*s
                    frequency = zeros(nALEtopics, 1);
                    % for each filterbox search most frequent cluster
                    for  xx = rx(x):rx(x) + 2*s
                        for  yy = ry(y):ry(y) + 2*s
                            for  zz = rz(z):rz(z) + 2*s
                                t = V(xx, yy, zz);
                                if t > 0
                                    frequency(t) = frequency(t) + 1;
                                end
                            end
                        end
                    end
                    % smooth the center voxel in the filterbox
                    t = V(rx(x) + s, ry(y) + s, rz(z) + s); 
                    if t > 0 
                        tmax = 0;
                        imax = 0;
                        for t = 1:nALEtopics
                            if tmax < frequency(t)
                                tmax = frequency(t);
                                imax = t;
                            end
                        end
                        if V(rx(x) + s, ry(y) + s, rz(z) + s) > 0
                            VV(rx(x) + s, ry(y) + s, rz(z) + s) = imax;
                        else
                            VV(rx(x) + s, ry(y) + s, rz(z) + s) = 0;
                        end
                    end
                end
            end
        end
    end
    V = VV;
    Title = [Title ' Smoothed(' int2str(s) ')'];
    end   

    % workaround for flatmap formatting issues
    VV = V;
    V(:,:,:,1) = VV(:,:,:);  
    % workaround for nifti formatting issues
    info.ImageSize = [77 96 79];
    info.PixelDimensions = [2 2 2];
    if contains(Title, 'MDTB'); 
        info.ImageSize = [141 95 87]; 
        info.PixelDimensions = [1 1 1];
        info.MultiplicativeScaling = 1;
    end
    info.Description = 'Modified using winner-take-all'; 
    % nifti
    niftiwrite(V, [ALEfolder 'TEMP.nii'], info)
    copyfile ([ALEfolder 'TEMP.nii'], [ALEfolder GraphName ALEfile ' ' Title '.nii'])  
    % flatmap
    figure();
    ff = [ALEfolder 'TEMP.nii'];
    if contains(Title, 'MDTB')
        Data = suit_map2surf(ff, 'space', 'SUIT', 'stats', @mode);
    else
        Data = suit_map2surf(ff, 'space', 'SPM', 'stats', @mode);
    end
    suit_plotflatmap(Data, 'type', 'label', 'cmap', Colormap, 'bordersize', 4);
    title([GraphName ALEfile ' ' Title ]);
    colormap(Colors);
    colorbar('Ticks', Colors2Map(:, 1)/size(Topics, 1), 'TickLabels', Topics);
    exportgraphics(gcf, [ALEfolder GraphName ' Flatmap ' ALEfile ' ' Title '.png'])
end

delete([ALEfolder 'TEMP.nii']);
delete([ALEfolder 'TEMPcolormap.txt']);
end
  


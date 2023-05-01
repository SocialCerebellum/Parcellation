function dwarfs_val_winner
% Validation of winner-take-all parcellation
% For the two ALE analyses of each Train and Test data part, the
% winner-take-all clustering is computed and compared using Pearson and 
% Rand coefficients. This procedure is repeated for each cross-validation.

%
% ENTER TITLE, CLUSTERS AND THEIR LABELS HERE BELOW
%
Title =    'CE+WB Complete (#98) on Z (euclidean)';							
Clusters = [1	5  15	0 0;   2  3  18 19 0;  6   7  16  20 13; 11  0	0	0 0;   ...
            4   0	0	0 0;  10  0	 0	 0 0;  8	0	0	0 0; 14	21	0	0 0;  9	0	0	0 0; ...
            22 12  17   0 0];	
Labels   = {'empathy',      'conflict',     'decision',       'face',    ...
            'emotion',      'objects',      'encoding',       'attention',     'word', ...
            'sensorimotor'};	
%
%Title =    'CE+WB Complete (#98) on Z (euclidean) 7-networks';							
%Clusters = [1	5  15	0 0 0 0;   2  3 18 19 0 0 0;  6 7 16 20 13 10 8;   ...
%            4  11	0	0 0 0 0;  14 21	 0  0 0 0 0;  9	0  0  0  0  0 0; ...
%            22 12  17   0 0 0 0];	
%Labels   = {'empathy',         'conflict',         'decision', ...
%            'emotion',         'attention',        'word', ...
%            'sensorimotor'};	

% create inputbox
answer = inputdlg({'Unique Topics filename string (before .nii):', ...
    'Topics folder:', ...
    'Colormap Filename (in ALE folder):', ...
    'Smoothing filter (1, 2, ... or empty):', ...
    'Title:' ...
    'Percentage for Testing (%)', ...
    'Start of Replications', 'End of Replications', ...
    'Whithing steps for non-overlapping voxels across boundaries (eg., 5)', ...
    'Whithening Threshold (eg., from 30 % non-overlap) up to + 50% in steps of 5%', ...
    'Figures/Flatmaps of each Replication (y or n)', ...
    'Full *.nii Filename with same Clustering and Smoothing'}, ...
    'Input', 1, ...
    {'_Z', 'G:\CerebellumDwarfs\4. CROSSVAL', 'Buckner-likeColors.txt', '2', Title, ...
    '50', '1', '20', '5', '30', 'n', 'Clusters on _Z CE+WB Complete (#98) on Z (euclidean) CE Data Rescaled  Smoothed(2)'});
%    {'_Z', 'G:\CerebellumDwarfs\4. CROSSVAL', 'Buckner-likeColors.txt', '2', Title, ...
%    '5', '1', '20', '3', '3', 'n', 'Clusters on _Z CE+WB Complete (#98) on Z (euclidean) CE Data Rescaled  Smoothed(2)'});
%    {'_Z', 'G:\CerebellumDwarfs\4. CROSSVAL', 'Buckner-likeColors.txt', '', Title, ...
%    '5', '1', '20', '1', '', 'n', 'Clusters on _Z CE+WB Complete (#98) on Z (euclidean) CE Data Rescaled '});


ALEfolder = [answer{2} '\'];
Colormap = answer{3};
Smoothing = answer{4};
Title = answer{5};
VALtestPc = str2num(answer{6});
VALrepsta = str2num(answer{7});
VALrepend = str2num(answer{8});
nBoundaryColors = str2num(answer{9});
Threshold = str2num(answer{10});
Flatmaps = answer{11};
if ~isempty(answer{12})
    ClusterFile = [answer{12} '.nii'];
else
    ClusterFile = '';
end
    info = niftiinfo([ALEfolder ClusterFile]);

% defaults
nClusters = size(Clusters, 1);
Colormap = [ALEfolder Colormap];
bColormap = Colormap; % for boundaries

XC = [];
% Write on Excel
xlsfilename = [ALEfolder 'Pearson & Rand Validation - ' Title '.xls'];
sheet = 1;
xlswrite(xlsfilename, {'replication'}, sheet, 'A1');


%=========
% Colormap
%=========

[CM(:, 1) CM(:, 2) CM(:, 3) CM(:, 4) CMlabels] = readvars(Colormap);
Colors = [CM(:,2) CM(:,3) CM(:,4)]/255;
Topics = char(CMlabels);

% Reordering of cluster colors for 'Buckner-likeColors.txt'
if contains(answer{3}, 'Buckner-likeColors.txt')
    CC = [];
    for i = 1:size(Topics, 1)
        for j = 1:size(Labels, 2)
            s1 = Topics(i, :);
            s2 = char(Labels(j));
            if contains(s1, s2)
                CC(j,:) = Colors(i,:); 
            end
        end
    end
    Colors = CC;
end
Topics = Labels';

% Add colors for whitening / transparancy
nColors = size(Colors, 1);
for i = 1:nBoundaryColors
    for j = 1:size(Colors, 1)
        Colors((nColors * i) + j, :) = Colors(j, :) * (i + 1);
    end
end

for i = 1:size(Colors, 1)
    Colors2Map(i, :) = i;
    Colors2Map(i, [2:4]) = Colors(i, [1:3]) * 255;
end

save([ALEfolder 'TEMPcolormap.txt'], 'Colors2Map', '-ascii');
Colormap = [ALEfolder 'TEMPcolormap.txt'];

%=======================
% START CROSS-VALIDATION
%=======================

for C1 = VALrepsta:VALrepend
    disp(['=========================']);
    disp(['== Cross-validation ' num2str(C1) ' ==']);
    disp(['=========================']);

    for T = 1:2
        if T == 1; VAL = 'train'; end
        if T == 2; VAL = 'test'; end
        
        ALEfilestring = ['Topic #' '*_C' VAL num2str(C1) '_' num2str(VALtestPc) '%' answer{1} '*.nii'];
        Title = answer{5};
        Title = [Title ' C' VAL num2str(C1)];

        % defaults
        nALEtopics = 0;
        ALEpath = [ALEfolder ALEfilestring];

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
        %           ~(contains(ALEfilenames(T1).name, '#0')) & ...
                if ~(contains(ALEfilenames(T1).name, 'ALE')) & ...
                   ~(contains(ALEfilenames(T1).name, '_pID')) & ...
                   ~(contains(ALEfilenames(T1).name, 'TEMP')) & ...
                   ~(contains(ALEfilenames(T1).name, 'Topics on')) & ...
                   ~(contains(ALEfilenames(T1).name, 'NoTopic')) & ...
                   ~(contains(ALEfilenames(T1).name, '#99')) & ...
                   ~(contains(ALEfilenames(T1).name, '#0')) & ...
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
        disp(['Cross-validation = ' num2str(C1)]);
        if nALEtopics == 0
            VALrepend = C1 - 1;
            break
        end

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

        ALEfile = [' on ' answer{1}(1:size(answer{1}, 2))];

        % name of graph
        GraphName = 'Clusters'; 
        % rewrite from linear X to volume matrix V
        V = niftiread(info);
        i = 0;
        for x = 1:size(V, 1)
            for y = 1:size(V, 2)
                for z = 1:size(V, 3)
                    i = i + 1;
                    V(x, y, z) = single(XmaxCluster(i)); 
                end
            end
        end

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

        % SUIT resampling
        % workaround for flatmap formatting issues
        VV = V;
        V(:,:,:,1) = VV(:,:,:);  
        % workaround for nifti formatting issues
        info.ImageSize = [77 96 79];
        %info.ImageSize = [80 96 70];
        info.PixelDimensions = [2 2 2];
        info.Description = 'Modified using winner-take-all'; 
        % nifti
        niftiwrite(V, [ALEfolder 'TEMP.nii'], info)
        copyfile ([ALEfolder 'TEMP.nii'], [ALEfolder GraphName ALEfile ' ' Title '.nii'])  
        % flatmap
        spm_suit
        ff = [ALEfolder 'TEMP.nii'];
        Data = suit_map2surf(ff, 'space', 'SPM', 'stats', @mode);
        % output graphics
        if Flatmaps ~= 'n'
            figure();
            suit_plotflatmap(Data, 'type', 'label', 'cmap', Colormap, 'bordersize', 4);
            title([GraphName ALEfile ' ' Title ]);
            colormap(Colors(1:nColors,:));
            colorbar('Ticks', Colors2Map(1:nColors, 1)/size(Topics, 1), 'TickLabels', Topics);
            exportgraphics(gcf, [ALEfolder GraphName ' Flatmap ' ALEfile ' ' Title '.png'])
        end
        % for visualization of overlapping voxels
        if T == 1; XC((C1 * 2 - 1),:) = XmaxCluster; end
        if T == 2; XC((C1 * 2),:) = XmaxCluster; end

        % original voxels
        if T == 1; Xtrain = XmaxCluster; end
        if T == 2; Xtest = XmaxCluster; end
        % use standardized data from SUIT with equal vector length of 28935 voxels
        for V1 = 1:size(Data)
            if isnan(Data(V1)); Data(V1) = 0; end
        end
        if T == 1; X1 = zeros(size(Data)); end 
        if T == 2; X2 = zeros(size(Data)); end 
        if T == 1; X1 = single(Data); end 
        if T == 2; X2 = single(Data); end 

    end % of T

    xlswrite(xlsfilename, {'replication'}, sheet, 'A1');
    xlswrite(xlsfilename, {'Orig: Pearson', 'p', 'Rand', 'Adj Rand'}, sheet, 'B1');
    xlswrite(xlsfilename, {'SUIT: Pearson', 'p', 'Rand', 'Adj Rand'}, sheet, 'AB1');

    row = C1 + 1;
    xlswrite(xlsfilename, C1, sheet, ['A' num2str(row)]);
    
    % original voxels
    XXtrain = Xtrain;
    XXtest = Xtest;
    for V = 1:size(XmaxCluster, 2)
        if Xtrain(V) == 0; XXtrain(V) = NaN; end
        if Xtest(V) == 0; XXtest(V) = NaN; end
    end
    [R,P] = corrcoef(XXtrain', XXtest', 'Rows','pairwise');
    xlswrite(xlsfilename, R(1, 2), sheet, ['B' num2str(row)]);
    xlswrite(xlsfilename, P(1, 2), sheet, ['C' num2str(row)]);
    
    V2 = 0;
    XXtrain = [];
    XXtest = [];
    for V = 1:size(XmaxCluster, 2)
        if Xtrain(V) ~= 0 & Xtest(V) ~= 0 
            V2 = V2 + 1;
            XXtrain(V2) = Xtrain(V);
            XXtest(V2) = Xtest(V);
        end
    end
    ri = rand_index(XXtrain, XXtest);
    ari = rand_index(XXtrain, XXtest, 'adjusted');
    xlswrite(xlsfilename, ri, sheet, ['D' num2str(row)]);
    xlswrite(xlsfilename, ari, sheet, ['E' num2str(row)]);
    
    % Rand per cluster
    xlswrite(xlsfilename, {'Cluster:'}, sheet, 'F1');
    for C = 1:nClusters
        XX1 = XXtrain;
        XX2 = XXtest;
        for V1 = 1:size(X1)
        % 9 (= 'other' category)
        if XXtrain(V1) ~= C; XX1(V1) = 9; end 
        if XXtest(V1) ~= C; XX2(V1) = 9; end
        end
    ri = rand_index(XX1, XX2);
    ari = rand_index(XX1, XX2, 'adjusted');
    xlswrite(xlsfilename, Labels(C), sheet, [char(64 + 6 + C) '1']);
    xlswrite(xlsfilename, Labels(C), sheet, [char(64 + 16 + C) '1']);
    xlswrite(xlsfilename, ri,  sheet, [char(64 + 6 + C) num2str(row)]);
    xlswrite(xlsfilename, ari, sheet, [char(64 + 16 + C) num2str(row)]);
    end

    % SUIT resampling
    [R,P] = corrcoef(X1, X2);
    xlswrite(xlsfilename, R(1, 2), sheet, ['AB' num2str(row)]);
    xlswrite(xlsfilename, P(1, 2), sheet, ['AC' num2str(row)]);
    ri = rand_index(X1, X2);
    ari = rand_index(X1, X2, 'adjusted');
    xlswrite(xlsfilename, ri, sheet, ['AD' num2str(row)]);
    xlswrite(xlsfilename, ari, sheet, ['AE' num2str(row)]);
    
    % Rand per cluster
    xlswrite(xlsfilename, {'Cluster:'}, sheet, 'AF1');
    for C = 1:nClusters
        XX1 = X1;
        XX2 = X2;
        for V1 = 1:size(X1)
        % 9 (= 'other' category)
        if X1(V1) ~= C; XX1(V1) = 9; end 
        if X2(V1) ~= C; XX2(V1) = 9; end
        end
    ri = rand_index(XX1, XX2);
    ari = rand_index(XX1, XX2, 'adjusted');
    xlswrite(xlsfilename, Labels(C), sheet, ['A' char(64 + 6 + C) '1']);
    xlswrite(xlsfilename, Labels(C), sheet, ['A' char(64 + 16 + C) '1']);
    xlswrite(xlsfilename, ri,  sheet, ['A' char(64 + 6 + C) num2str(row)]);
    xlswrite(xlsfilename, ari, sheet, ['A' char(64 + 16 + C) num2str(row)]);
    end

end % end of C1

%=====================
% END CROSS-VALIDATION
%=====================

Title = [answer{5}];

% Full original cluster nii file (10 clusters)
if ~isempty(ClusterFile)
    XF = [];
    info = niftiinfo([ALEfolder ClusterFile]);
    V = niftiread(info);
    % write in linear (volume) matrix FX
    V1 = 0;
    for x = 1:size(V, 1)
        for y = 1:size(V, 2)
            for z = 1:size(V, 3)
                V1 = V1 + 1;
                XF(V1) = V(x, y, z);
            end
        end
    end     
end

% Convert cluster coding to 7 networks
if nClusters == 7
    Conversions =  [1	1;	2	2;	3	3;	4	4;	5	4;	...
                    6	3;	7	3;	8	5;	9	6;	10	7];	
     for V1 = 1:size(XF, 2)
        if XF(V1) > 0 
            for c = 1:size(Conversions)
                if XF(V1) == Conversions(c, 1)
                   XF(V1) = Conversions(c, 2);
                   break
                end
            end
        end
    end
end

% Compute misplaced voxels
if ~isempty(ClusterFile)  
    % number of misplaced voxels
    for T2 = 1:size(XC, 2)
        Xmspl(T2) = 0;
        for T3 = 1:size(XC, 1)
            if XF(T2) ~= XC(T3, T2)
                Xmspl(T2) = Xmspl(T2) + 1;
            end 
        end
    end
end

% Redraw colors according to misplacements
for TH1 = Threshold:5:(Threshold + 50)
    if TH1 > 100
        break
    end
    TH = TH1 / 100;
    XXF = XF;
    for T2 = 1:size(XC, 2)
        if XF(T2) > 0 & Xmspl(T2) > 0
            %XXF(T2) 
            %Xmspl(T2)
            if Xmspl(T2) / size(XC, 1) > TH
                % divide misplaced voxels by number of cross-validations to obtain % misplacements
                % subtract the % threshold (eg 0.30)
                % divide that subtraction by (1 - threshold; 0.70) to resize to full 100%
                % multiply all of this by number of color boundary steps 
                bb = round((Xmspl(T2) / size(XC, 1) - TH) / (1 - TH) * nBoundaryColors + 1);
                if bb > nBoundaryColors; bb = nBoundaryColors; end
                XXF(T2) = XF(T2) + (bb * nColors);
                %XXF(T2)
                %XXF(T2)
            end
        end  
    end

    % Create Flatmap
    GraphName = ['Validation threshold ' num2str(TH *100) '%']; 
    % rewrite from linear X to volume matrix V
    V = niftiread(info);
    V1 = 0;
    for x = 1:size(V, 1)
        for y = 1:size(V, 2)
            for z = 1:size(V, 3)
                V1 = V1 + 1;
                V(x, y, z) = single(XXF(V1)); 
            end
        end
    end 

    % workaround for flatmap formatting issues
    VV = V;
    V(:,:,:,1) = VV(:,:,:);  
    % workaround for nifti formatting issues
    info.ImageSize = [77 96 79];
    %info.ImageSize = [80 96 70];
    info.PixelDimensions = [2 2 2];
    info.Description = 'Modified Validation Boundaries'; 
    % nifti
    niftiwrite(V, [ALEfolder 'TEMP.nii'], info)
    %copyfile ([ALEfolder 'TEMP.nii'], [ALEfolder GraphName ALEfile ' ' Title '.nii'])  
    % flatmap
    spm_suit
    figure();
    ff = [ALEfolder 'TEMP.nii'];
    Data = suit_map2surf(ff, 'space', 'SPM', 'stats', @mode);
    suit_plotflatmap(Data, 'type', 'label', 'cmap', Colormap, 'bordersize', 4);
    title([GraphName ALEfile ' ' Title]);
    colormap(Colors(1:nColors,:));
    colorbar('Ticks', Colors2Map(1:nColors, 1)/size(Topics, 1), 'TickLabels', Topics);
    exportgraphics(gcf, [ALEfolder GraphName ' Flatmap' ALEfile ' ' Title '.png'])
    end
end

function ri = rand_index(p1, p2, varargin)
%RAND_INDEX Computes the rand index between two partitions.
%   RAND_INDEX(p1, p2) computes the rand index between partitions p1 and
%   p2. Both p1 and p2 must be specified as N-by-1 or 1-by-N vectors in
%   which each elements is an integer indicating which cluster the point
%   belongs to.
%
%   RAND_INDEX(p1, p2, 'adjusted') computes the adjusted rand index
%   between partitions p1 and p2. The adjustment accounts for chance
%   correlation.

    % Parse the input and throw errors
    % Check inputs
    adj = 0;
    if nargin == 0
        error('Arguments must be supplied.');
    end
    if nargin == 1
        error('Two partitions must be supplied.');
    end
    if nargin > 3
        error('Too many input arguments');
    end
    if nargin == 3
        if strcmp(varargin{1}, 'adjusted')
            adj = 1;
        else
            error('%s is an unrecognized argument.', varargin{1});
        end
    end
    if length(p1)~=length(p2)
        error('Both partitions must contain the same number of points.');
    end

    % Check if arguments need to be flattened
    if length(p1)~=numel(p1)
        p1 = p1(:);
        warning('The first partition was flattened to a 1D vector.')
    end
    if length(p2)~=numel(p2)
        p2 = p2(:);
        warning('The second partition was flattened to a 1D vector.')
    end

    % Check for integers
    if isreal(p1) && all(rem(p1, 1)==0)
        % all is well
    else
        warning('The first partition contains non-integers, which may make the results meaningless. Attempting to continue calculations.');
    end

    if isreal(p2) && all(rem(p2, 1)==0)
        % all is well
    else
        warning('The second partition contains non-integers, which may make the results meaningless. Attempting to continue calculations.');
    end

	% Preliminary computations and cleansing of the partitions
    N = length(p1);
    [~, ~, p1] = unique(p1);
    N1 = max(p1);
    [~, ~, p2] = unique(p2);
    N2 = max(p2);

    % Create the matching matrix
    for i=1:1:N1
        for j=1:1:N2
            G1 = find(p1==i);
            G2 = find(p2==j);
            n(i,j) = length(intersect(G1,G2));
        end
    end

    % If required, calculate the basic rand index
    if adj==0
        ss = sum(sum(n.^2));
        ss1 = sum(sum(n,1).^2);
        ss2 = sum(sum(n,2).^2);
        ri = (nchoosek2(N,2) + ss - 0.5*ss1 - 0.5*ss2)/nchoosek2(N,2);
    end


    % Otherwise, calculate the adjusted rand index
    if adj==1
        ssm = 0;
        sm1 = 0;
        sm2 = 0;
        for i=1:1:N1
            for j=1:1:N2
                ssm = ssm + nchoosek2(n(i,j),2);
            end
        end
        temp = sum(n,2);
        for i=1:1:N1
            sm1 = sm1 + nchoosek2(temp(i),2);
        end
        temp = sum(n,1);
        for i=1:1:N2
            sm2 = sm2 + nchoosek2(temp(i),2);
        end
        NN = ssm - sm1*sm2/nchoosek2(N,2);
        DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);

        % Special case to handle perfect partitions
        if (NN == 0 && DD==0)
            ri = 1;
        else
            ri = NN/DD;
        end
    end


    % Special definition of n choose k
    function c = nchoosek2(a,b)
        if a>1
            c = nchoosek(a,b);
        else
            c = 0;
        end
    end
end  
  


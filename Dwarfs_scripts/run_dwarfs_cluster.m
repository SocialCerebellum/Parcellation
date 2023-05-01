function dwarfs_cluster

% create inputbox
answer = inputdlg({'Unique ALE filename string (before .nii):', ...
    'ALE folder:', 'Include MDTB (y, n, o[nly])', ...
    'Cluster (kmeans, ward, complete, weighted, average, all)', ...
    'Cluster metric (euclidean [sq for kmeans], sqeuclidean, correlation [only kmeans])', ...
    'Cluster maximum (or leave empty)', ...
    'Create Flatmaps (y, n, o[nly])'}, ...
    'Input', 1, ...
    {'Set', 'G:\CerebellumDwarfs\MDTB', 'n', 'all', 'euclidean', '15', 'n'});
%    {'_Z', 'G:\CerebellumDwarfs\3. ALE_CE_WB', 'n', 'all', 'euclidean', '20', 'n'}); 
%    {'Set', 'G:\CerebellumDwarfs\MDTB', 'n', 'all', 'euclidean', '15', 'n'});
%    {'_Z', 'G:\CerebellumDwarfs\3. ALE', 'n', 'all', 'euclidean', '15', 'n'});
%    {'_Z', 'G:\CerebellumDwarfs\3. ALE_WB', 'n', 'all', 'euclidean', '15', 'n'});

ALEfilestring = ['*' answer{1} '*.nii'];
ALEfolder = [answer{2} '\'];
IncludeMDTB = answer{3};
ClusterMethod = answer{4};
ClusterMetric = answer{5};
ClusterMaximum = str2num(answer{6});  
Flatmaps = answer{7};

MDTBfolder = 'G:\CerebellumDwarfs\MDTB\';

nALEtopics = 0;
maxV_ALE = 0;
maxV_MDTB = 0;
minV_ALE = 0;
minV_MDTB = 0;
cluster = '0';

ALEpath = [ALEfolder ALEfilestring];

% activate SUIT toolbox
if Flatmaps ~= 'n'
    spm_suit
end

%========================
% Select and Read Volumes
%========================

% empty ppdirectories
ff = 0;
ffdirectories = {};
% for each filename
TT1 = 0;
% make a list of all ALE files:
if IncludeMDTB ~= 'o'
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
    disp(['Number of Topics = ' num2str(nALEtopics)])
end
    
if IncludeMDTB ~= 'n'
    MDTBfilenames = dir([MDTBfolder '*.nii']);
    % for each filename
    for T3 = 1:(size(MDTBfilenames, 1))
        try
            MDTBdirectories{T3}= MDTBfilenames(T3).name;
            ffdirectories{nALEtopics + T3} = MDTBfilenames(T3).name;
        continue
        end
    end
    ff = ff + T3;
end

% read volumes
for T2 = 1:ff
    T2s(T2) = T2;

    if T2 <= nALEtopics
        if IncludeMDTB ~= 'o'
            ALEfilename = ffdirectories{T2};
            V = niftiread([ALEfolder ALEfilename]);
            % for rescaling MDTB
            if minV_ALE > min(V, [], 'all')
                minV_ALE = min(V, [], 'all'); 
            end
            if maxV_ALE < max(V, [], 'all')
                maxV_ALE = max(V, [], 'all'); 
            end
        end
    else
        T3 = T2 - nALEtopics;
        MDTBfilename = MDTBdirectories{T3};
        V = niftiread([MDTBfolder MDTBfilename]);
        % for rescaling MDTB
        if minV_MDTB > min(V, [], 'all')
            minV_MDTB = min(V, [], 'all'); 
        end
        if maxV_MDTB < max(V, [], 'all')
            maxV_MDTB = max(V, [], 'all'); 
        end 
    end
     
    %write in linear (volume) matrix X
    i = 0;
    for x = 1:size(V, 1)
        for y = 1:size(V, 2)
            for z = 1:size(V, 3)
                i = i + 1;
            X(T2, i) = single(V(x, y, z));
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
        s = ALEfilename(1:size(ALEfilename, 2) - 4);
        title(s);
        exportgraphics(gcf, [ALEfolder 'Flatmap' s '.png'])
    end

end

if Flatmaps == 'o'
    return
end

%========================
% Reformat for Clustering
%========================

% ommit NaN columns from linear volume matrix
ii = 0;
for i = 1:size(X, 2)
    isum = 0;
    for T2 = 1:size(X, 1)
        isum = isum + X(T2, i);
    end
    if isum > 0 % voxels are not all empty
        ii = ii + 1;
        XX(1:size(X, 1), ii) = X(1:size(X, 1), i);
    end
end

% rescale MDTB
for T2  = nALEtopics + 1:ff 
    for i = 1:size(XX, 2)
        % intercept 
%OFF%        XX(T2, i) = XX(T2, i) - single(minV_MDTB + minV_ALE);
        % range
%OFF%        XX(T2, i) = XX(T2, i) / single(maxV_MDTB - minV_MDTB) * single(maxV_ALE - minV_ALE);
    end
end

% ==========
% clustering
% ==========

Metric = ClusterMetric;
answer{1} = [answer{1} ' ' Metric];
if isempty(ClusterMaximum) 
    ClusterMaximum = nALEtopics;
end

if IncludeMDTB == 'y' 
    nALEtopics = nALEtopics + T3;
    answer{1} = [answer{1} ' (+ MDTB)'];
elseif IncludeMDTB == 'o'
    nALEtopics = nALEtopics + T3;
    answer{1} = [answer{1} ' (only MDTB)'];
end

xlsfilename = [ALEfolder 'Clusters on ' answer{1} '.xls'];

if contains(ClusterMethod,'kmeans') | contains(ClusterMethod,'all')
    cluster = 'kmeans';
    disp('Running clustering kmeans')
    sheet = 1;
    kmeansMetric = ClusterMetric;
    if contains(ClusterMetric, 'euclidean')
        kmeansMetric = 'sqeuclidean';
    end
    for T2 = 1:min(nALEtopics, ClusterMaximum)
        disp (['   Clusters = ' num2str(T2)])
        opts = [];
        idx = kmeans(XX,T2,'Replicates',10,'Distance', kmeansMetric, 'Options',opts);
        idxs(:, T2) = idx(:);
    end
    i = min(nALEtopics, ClusterMaximum);
    eva = evalclusters(XX,'kmeans','CalinskiHarabasz','KList',1:i);
    figure()
    plot(eva)
    title([cluster ' on ' answer{1}])
    exportgraphics(gcf,[ALEfolder cluster ' on ' answer{1} '.png' ])
    xlswrite(xlsfilename, {['Method: '] cluster}, sheet, 'A1');
    xlswrite(xlsfilename, {'Number of Clusters >>>>'}, sheet, 'C1');
    xlswrite(xlsfilename, T2s(1:min(nALEtopics, ClusterMaximum)), 1, 'F1');
    xlswrite(xlsfilename, T2s', sheet, 'A3');
    xlswrite(xlsfilename, ffdirectories', sheet, 'B3');
    xlswrite(xlsfilename, idxs, sheet, 'F3');
end
if ~contains(ClusterMethod,'kmeans')
for i = 1:4
    if i == 1
        if contains(ClusterMethod,'ward') | contains(ClusterMethod,'all')
            cluster = 'ward';
        end
    elseif i == 2
        if contains(ClusterMethod,'complete') | contains(ClusterMethod,'all')
            cluster = 'complete';
        end  
    elseif i == 3
        if contains(ClusterMethod,'weighted') | contains(ClusterMethod,'all')
            cluster = 'weighted';
        end
    elseif i == 4
        if contains(ClusterMethod,'average') | contains(ClusterMethod,'all')
            cluster = 'average';
        end
    end 
    if ~(contains(cluster,'0'))
        disp(['Running clustering ' cluster])
        sheet = i + 1;
        end
        if contains(ClusterMetric, 'sqeuclidean')
            Metric = 'squaredeuclidean';
        end
        for i = 1:min(nALEtopics, ClusterMaximum);
            R = clusterdata(XX, 'Linkage', cluster, 'Maxclust', i, 'Distance', Metric);
            RR(1:size(R, 1), i) = R;
        end;
        tree = linkage(XX, cluster, Metric);
        figure()
        [H,T,outperm] = dendrogram(tree, 0, 'Orientation', 'left');
        title([cluster ' on ' answer{1}])
        exportgraphics(gcf,[ALEfolder cluster ' on ' answer{1} '.png'])
        xlswrite(xlsfilename, {['Method: '] cluster}, sheet, 'A1');
        xlswrite(xlsfilename, T2s', sheet, 'A3');
        xlswrite(xlsfilename, ffdirectories', sheet, 'B3');
        xlswrite(xlsfilename, {'Number of Clusters >>>>'}, sheet, 'I1');
        xlswrite(xlsfilename, R, sheet, 'E3');
        xlswrite(xlsfilename, {'Tree'}, sheet, 'G1');
        xlswrite(xlsfilename, outperm', sheet, 'G3');
        for i = 1:nALEtopics 
            HH(i) = ffdirectories(outperm(i));
            for ii = 1:ClusterMaximum
                LL(i, ii) = RR(outperm(i),ii);
            end
        end   
        xlswrite(xlsfilename, {'Tree'}, sheet, 'H1');
        xlswrite(xlsfilename, HH', sheet, 'H3');
        xlswrite(xlsfilename, T2s(1:min(nALEtopics, ClusterMaximum)), sheet, 'L1');
        xlswrite(xlsfilename, LL, sheet, 'L3');            
        cluster = '0';
    end
end

end
    


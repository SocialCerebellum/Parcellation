function dwarfs_rand
% Rand index to compare the parcellations 

% Clusters can be converted according to Buckner-like clusters 
% ENTER CLUSTER CONVERSIONS FURTHER BELOW IN CODE

% ENTER PARCELLATION FILES HERE BELOW

Parcellations = [{'Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii',...	
                'Ji_10Networks.nii',...	
                'MDTB_10Regions.nii',...	
                'Clusters on Set Ward 61 on MDTB (euclidean) Rescaled  Smoothed(2).nii',...
                'Clusters on _Z CE+WB Complete (#98) on Z (euclidean) CE Data Rescaled  Smoothed(2).nii'}];

% create inputbox
answer = inputdlg({'Folder:', 'ColorFile', '(Flat)maps after cluster conversion (y or n)'}, ...
    'Input', 1, ...
    {'G:\CerebellumDwarfs', 'Buckner2011_7+1Networks_Color.txt', 'n'}); 
    
Folder = [answer{1} '\'];

% loop through all parcellations (lower portion of matrix)
for P1 = 1:size(Parcellations, 2) - 1
    for P2 = P1 + 1:size(Parcellations, 2)

    NiftiFile1 = char(Parcellations(P1))
    NiftiFile2 = char(Parcellations(P2))
    ColorFile = answer{2};
    ColorFile = [Folder ColorFile];
    Flatmaps = answer{3};

    %=========
    % Colormap
    %=========

    [CM(:, 1) CM(:, 2) CM(:, 3) CM(:, 4) CMlabels] = readvars(ColorFile);
    Colors = [CM(:,2) CM(:,3) CM(:,4)]/255;
    Topics = char(CMlabels);
    for i = 1:size(Colors, 1)
        Colors2Map(i, :) = i;
        Colors2Map(i, [2:4]) = Colors(i, [1:3]) * 255;
    end
    save([Folder 'TEMPcolormap.txt'], 'Colors2Map', '-ascii');
    ColorFile = [Folder 'TEMPcolormap.txt'];


    %================================
    % Convert clusters & Read Volumes
    %================================

    % loop through each of the two parcellations that are compared
    for T = 1:2
        % read volumes
        X = [];
        if T == 1; NiftiFile = NiftiFile1; end
        if T == 2; NiftiFile = NiftiFile2; end

        % defaults
        Parcellation = NiftiFile;
        Titlemap = NiftiFile(1:length(NiftiFile)-4);
        NiftiFile = [Folder NiftiFile];
        Conversions = [];

        %===============================
        % ENTER CLUSTER CONVERSIONS HERE
        %===============================

        if contains(Parcellation, 'Ji_10Networks.nii')	
            Conversions =   [1	1;	...
                            2	1;	...
                            3	2;	...
                            4	4;	...
                            5	3;	...
                            6	8;	...
                            7	6;	...
                            8	9;	...
                            9	7;	...
                            10	9];	
        elseif contains(Parcellation, 'MDTB_10Regions.nii')	
            Conversions =   [1	2
                            2	2;	...
                            3	1;	...
                            4	9;	...
                            5	4;	...
                            6	4;	...
                            7	8;	...
                            8	8;	...
                            9	8;	...
                            10	9];	
        elseif contains(Parcellation, 'Clusters on Set Ward 61 on MDTB (euclidean) Rescaled  Smoothed(2).nii')	
            Conversions =   [1	7;	...
                            2	8;	...
                            3	6;	...
                            4	7;	...
                            5	2;	...
                            6	3;	...
                            7	6;	...
                            8	9;	...
                            9	7;	...
                            10	4;	...
                            11	5;	...
                            12	6;	...
                            13	9;	...
                            14	6;	...
                            15	4];	
        elseif contains(Parcellation, 'Clusters on _Z CE+WB Complete (#98) on Z (euclidean) CE Data Rescaled  Smoothed(2).nii')	
            Conversions =   [1	7;	...
                            2	4;	...
                            3	6;	...
                            4	5;	...
                            5	5;	...
                            6	6;	...
                            7	6;	...
                            8	3;	...
                            9	8;	...
                            10	2];	
        end

        info = niftiinfo(NiftiFile);
        V = niftiread(info);
        % write in linear (volume) matrix X
        V1 = 0;
        for x = 1:size(V, 1)
            for y = 1:size(V, 2)
                for z = 1:size(V, 3)
                    V1 = V1 + 1;
                    X(T, V1) = V(x, y, z);
                end
            end
        end

        % convert cluster coding
        for V1 = 1:size(X, 2)
            if isnan(X(T, V1)); X(T, V1) = 0; end
            if X(T, V1) > 0 & size(Conversions) > 0
                for c = 1:size(Conversions)
                    if X(T, V1) == Conversions(c, 1)
                       X(T, V1) = Conversions(c, 2);
                       break
                    end
                end
            end
        end

        % rewrite from linear X to volume matrix V
        V = niftiread(info);
        i = 0;
        for x = 1:size(V, 1)
            for y = 1:size(V, 2)
                for z = 1:size(V, 3)
                    i = i + 1;
                    V(x, y, z) = X(T, i); 
                end
            end
        end
        %Vtest= V(20:40, 20:35, 40:50)  


        %========
        % flatmap
        %========

        % workaround for flatmap formatting issues
        VV = V;
        V(:,:,:,1) = VV(:,:,:);  
        % workaround for nifti formatting issues
        info.ImageSize = [77 96 79];
        info.PixelDimensions = [2 2 2];

        if contains(Parcellation, 'Buckner2011_7Networks_MNI152')
            info.ImageSize = [256 256 256]; 
            info.PixelDimensions = [1 1 1];
            info.MultiplicativeScaling = 1;
        end
        if contains(Parcellation, 'Ji_10Networks')
            info.ImageSize = [141 95 87]; 
            info.PixelDimensions = [1 1 1];
            info.MultiplicativeScaling = 1;
        end
        if contains(Parcellation, 'MDTB_10Regions.nii')
            info.ImageSize = [71 48 44]; 
            info.PixelDimensions = [2 2 2];
            info.MultiplicativeScaling = 1;
        end
        if contains(Parcellation, 'Clusters on Set Ward 61 on MDTB'); 
            info.ImageSize = [141 95 87]; 
            info.PixelDimensions = [1 1 1];
            info.MultiplicativeScaling = 1;
        end
        info.Description = 'Modified using cluster conversion'; 

        % nifti
        niftiwrite(V, [Folder 'TEMP.nii'], info)
        copyfile ([Folder 'TEMP.nii'], [Folder 'Cluster conversion ' Parcellation])  
        spm_suit
        ff = [Folder 'TEMP.nii'];
        Data = suit_map2surf(ff, 'space', 'SPM', 'stats', @mode);

        %=========
        % plotting
        %=========
        if Flatmaps ~= 'n' & P1 == 1
            figure();
            suit_plotflatmap(Data, 'type', 'label', 'cmap', ColorFile, 'bordersize', 4);
            title(Titlemap);
            colormap(Colors);
            colorbar('Ticks', Colors2Map(:, 1)/size(Topics, 1), 'TickLabels', Topics);
            exportgraphics(gcf, [Folder 'Cluster conversion Flatmap ' Titlemap '.png'])
        end

        % use standardized data from SUIT with equal vector length of 28935 voxels
        for V1 = 1:size(Data)
            if isnan(Data(V1)); Data(V1) = 0; end
        end
        if T == 1; X1 = zeros(size(Data)); end 
        if T == 2; X2 = zeros(size(Data)); end 
        if T == 1; X1 = single(Data); end 
        if T == 2; X2 = single(Data); end 

    end % of T

    % Write on Excel
    xlsfilename = [Folder 'Adjusted Rand index.xls'];
    sheet = 1;

    ari = rand_index(X1, X2, 'adjusted')
    xlswrite(xlsfilename, ari, sheet, [char(P1 + 64 + 1) num2str(P2 + 1)]);

    ri = rand_index(X1, X2)
    xlswrite(xlsfilename, ri, sheet, [char(P1 + 64 + 7) num2str(P2 + 1)]);

    % Rand per cluster
    for C = 1:8
        XX1 = X1;
        XX2 = X2;
        for V1 = 1:size(X1)
        % assymetric convergence when using NaN instead of 9 (= 'other' category)
        if X1(V1) ~= C; XX1(V1) = 9; end 
        if X2(V1) ~= C; XX2(V1) = 9; end
        end
    ari = rand_index(XX1, XX2, 'adjusted')
    xlswrite(xlsfilename, ari, sheet, [char(P1 + 64 + 1) num2str(P2 + 1 + C * 7)]);

    ri = rand_index(XX1, XX2)
    xlswrite(xlsfilename, ri, sheet, [char(P1 + 64 + 7) num2str(P2 + 1 + C * 7)]);
    end


    end % of P2
end % of P1

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


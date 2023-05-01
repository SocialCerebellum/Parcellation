function dwarfs_recolor

% create inputbox
answer = inputdlg({'folder:', 'nifti file', ...
    'Colormap Filename (in folder):'}, ...
    'Input', 1, ...
    {'G:\CerebellumDwarfs', 'MDTB_10Regions.nii', ...
     'MDTB-10Regions-BucknerColor-Labeled.txt'});
%    {'G:\CerebellumDwarfs', 'Ji_10Networks.nii', ... 
%     'Ji-10Regions-BucknerColor-Labeled.txt'});
%    {'G:\CerebellumDwarfs', 'Dwarfs-10Clusters(2).nii', ...
%     'Dwarfs-10Regions-BucknerColor-Labeled.txt'});

Folder = [answer{1} '\'];
NiftiFile = answer{2};
ColorFile = answer{3};
Titlemap = ColorFile(1:length(ColorFile)-4);

% defaults
NiftiFile = [Folder NiftiFile];
ColorFile = [Folder ColorFile];

%=========
% Colormap
%=========

[CM(:, 1) CM(:, 2) CM(:, 3) CM(:, 4) CMlabels] = readvars(ColorFile);
Colors = [CM(:,2) CM(:,3) CM(:,4)]/255;
%Colors = [CM(:,2) CM(:,3) CM(:,4)];
Topics = char(CMlabels);
for i = 1:size(Colors, 1)
    Colors2Map(i, :) = i;
    Colors2Map(i, [2:4]) = Colors(i, [1:3]) * 255;
    %Colors2Map(i, [2:4]) = Colors(i, [1:3]);
end
save([Folder 'TEMPcolormap.txt'], 'Colors2Map', '-ascii');
ColorFile = [Folder 'TEMPcolormap.txt'];

%========
% flatmap
%========

spm_suit
figure();
Data = suit_map2surf(NiftiFile, 'space', 'SPM', 'stats', @mode);
suit_plotflatmap(Data, 'type', 'label', 'cmap', ColorFile, 'bordersize', 4);
title(Titlemap);
colormap(Colors);
colorbar('Ticks', Colors2Map(:, 1)/size(Topics, 1), 'TickLabels', Topics);
exportgraphics(gcf, [Folder 'Flatmap ' Titlemap '.png'])
end

  


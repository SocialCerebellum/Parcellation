function dwarfs_cluster

% create inputbox
answer = inputdlg({'Unique ALE filename string (before .nii):', 'ALE folder:'}, 'Input', 1, ...
    {'pID1.0E-4_ALE', 'G:\CerebellumDwarfs\3. ALE'});
ALEfilestring = ['*' answer{1} '*.nii'];
ALEfolder = [answer{2} '\'];

ALEpath = [ALEfolder ALEfilestring];
cd(ALEfolder)

% activate SUIT toolbox
spm_suit

%ALE_Filename = 'Topic #1 = 1 - trait_pID1.0E-4_ALE.nii';
%ALE_folder = 'G:\CerebellumDwarfs\3. ALE';
%ALE_folder = [ALE_folder '\'];

%Make a list of all ALE files:
ALEfilenames = dir(ALEpath);
ff = (size(ALEfilenames, 1));
    % empty ppdirectories
    ffdirectories = {};
    % for each filename
    for S1 = 1:ff
        try
            ffdirectories{S1}= ALEfilenames(S1).name;
        continue
        end
    end

ALEfile = [];
for S2 = 1:ff  
   
    ALEfilename = ffdirectories{S2};
%    disp(['Writing thresholded flatmap of ' ALEfolder ALEfilename]) 

    % needs activation of SUIT
    reset(gcf);

    %SPM_HOME = 'C:\SPM\spm12';
    %Parc = suit_map2surf([SPM_HOME '\toolbox\suit\atlasesSUIT\Buckner_7Networks.nii'],'stats',@mode);
    %suit_plotflatmap(Parc,'type','label','cmap',[SPM_HOME '\toolbox\suit\atlasesSUIT\Buckner_7Networks_Color.txt'])
    %hold on;

%    Data = suit_map2surf([ALEfolder ALEfilename], 'space','SPM','stats',@minORmax);
    %%suit_plotflatmap(Data,'underlay',[SPM_HOME '\toolbox\suit\atlasesSUIT\Buckner_7Networks.nii'],'cscale',[],'cmap',hot);
%    suit_plotflatmap(Data,'cscale',[],'cmap',hot);
%    title(ALEfilename);
    hold off;

    %exportgraphics(gcf,'test.png')
    %exportgraphics(gcf,ALEfilename)

    X = niftiread([ALEfolder ALEfilename]);
    
    %X = [randn(10,1)*0.75+ones(10,1)];
    %Y = [randn(10,1)*0.5-ones(10,1)];
    
    ALEfile = [ALEfile, X];
    
%    fff = ALEfile(1:10, 1:10, 1:10)
    
%    sz = size(ALEfile)
    %disp(ALEfile(1:10,1:S2))
    %disp(ALEfile((sz - 10):sz,1:S2))
end

rng('default');  % For reproducibility
X = [gallery('uniformdata',[10 3],12); ...
    gallery('uniformdata',[10 3],13)+1.2; ...
    gallery('uniformdata',[10 3],14)+2.5];
y = [ones(10,1);2*(ones(10,1));3*(ones(10,1))]; % Actual classes
Create a scatter plot of the data.
scatter3(X(:,1),X(:,2),X(:,3),100,y,'filled')
title('Randomly Generated Data in Three Clusters');

Find a maximum of three clusters in the data by specifying the value 3 for the cutoff input argument.
T1 = clusterdata(X,3);
Because the value of cutoff is greater than 2, clusterdata interprets cutoff as the maximum number of clusters.
Plot the data with the resulting cluster assignments.
scatter3(X(:,1),X(:,2),X(:,3),100,T1,'filled')
title('Result of Clustering');

Find a maximum of three clusters by specifying the value 3 for the 'MaxClust' name-value pair argument.
T2 = clusterdata(X,'Maxclust',3); 
Plot the data with the resulting cluster assignments.
scatter3(X(:,1),X(:,2),X(:,3),100,T2,'filled')
title('Result of Clustering');

return

rng default; % For reproducibility
X = [randn(100,100,2)*0.75+ones(100,100,2); randn(100,100,2)*0.5-ones(100,100,2)];

figure;
plot(X(:,1),X(:,2),'.');
title 'Randomly Generated Data';

opts = statset('Display','final');
%[idx,C] = kmeans(X,2,'Distance','cityblock','Replicates',5,'Options',opts);
[idx,C] = kmeans(X,2,'Replicates',5);

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids','Location','NW')
title 'Cluster Assignments and Centroids'
hold off

[silh,h] = silhouette(X,idx,'cityblock');
xlabel('Silhouette Value')
ylabel('Cluster')

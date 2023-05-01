# Parcellation
Parcellation of the cerebellum based on NeuroSynth database

Files related to the distinct step in the analysis


Selection of Cerebellar Coordinates

1.	Accessing NeuroSynth Database on 12 November 2020:
saved in Excel file database_select
2.	Selection of a cerebellum mask (ROI) by anatomical boundaries
•	Within a box defined by the outer x-y-z coordinates -63 > x < 63; -28 > y >-103; and z < 11
(Excel database_CE_ROI.xlsx in columns “box_x, box_y, box_z”)
•	Below the z coordinates surrounding the superior part of the cerebellum determined by two planes defined by x-y-z coordinates: 0 -50 12 and ±60 -50 -30, using the formulas: z < (-30 -12) / (60 -0) * ABS(x) +12; and z < -0.70 * ABS(x) +12
(Excel database_CE_ROI.xlsx in column “0 -50 12 / 60 -50 -30 plane”)
•	Below the z coordinates surrounding the declining posterior part of the cerebellum determined by x-y-z coordinates: 0 -50 12 and 0 -100 -15
(Excel database_CE_ROI.xlsx in column “0 -50 12 / 0 -100 -15 plane”)
3.	Selection using Brain Atlases:
•	Using the SPM Anatomy box 
Excel database_CE_ROI.xlsx in columns following “SPM ANATOMY”
(resulting in n = 26601 coordinates; n = 5726 articles
•	Using the Talairach Client toolbox
Excel database_CE_ROI.xlsx in columns following “TALAIRACH CLIENT”	
resulting in 29470 coordinates; 6177 articles
4.	Selection of only MNI coordinates:
Excel database_CE_pubs_tabs_MNI to be validated by the Dwarfs	
resulting in 4716 articles with MNI coordinates

Selection of Functional NeuroSynth Topics and Terms 

1.	Extracting Topics from NeuroSynth website
	https://neurosynth.org/analyses/topics/v5-topics-50/
saved in Excel file v5-topics-51.xlsx
2.	Selecting top 10 terms for each NeuroSynth topic 
Excel v5-topics-51.xlsx columns following “Top 10 terms”
3.	Excluding non-functional terms 
Excel v5-topics-51.xlsx column “Selected” switched to “0”
4.	Selecting only unique terms among the remaining functional terms
Excel v5-topics-51.xlsx columns “Uniqueness” parameter as low as possible and close to 1 and the columns under the heading “Unique terms top 10”
5.	Final selection of representative functional terms 
computed for overlap in Excel database_pubs_overlap_000_term and for uniqueness in Excel database_pubs_centrality_000_term.xlsx

Quality screening of Selected Studies and Coordinate Tables

1.	Excluding studies or tables referring to analyses with (see method section)
2.	Including studies or tables referring to analyses with (see method section)
3.	Combining all databases from individual collaborators (i.e., dwarfs)
database_CE_pubs_tabs_MNI_NOterms_AllDwarfs
database_CE_pubs_tabs_MNI_NOterms_AllDwarfs_ValidityCheck
database_CE_pubs_tabs_MNI_NOterms_AllDwarfs_ValidityCheck_Sorted
database_CE_pubs_tabs_MNI_NOterms_AllDwarfs_Valid
database_CE_pubs_tabs_Valid for all validated studies and tables
4.	Adding number of Participants:
database_CE_pubs_tabs_MNI_NOterms_Valid with 50325 participants
5.	Adding coordinates
Excel database_CE_valid validated for all 1873 studies and 12485 coordinates taken from the whole CE selected database

Functional Parcellation

1.	Preparing input for ALE Analysis of functional topics
MatLab: run_dwarfs_ALEinput
2.	ALE Analysis of individual topics
GingerAle run manually using the graphic interface
3.	Clustering of functional topics
MatLab: run_dwarfs_cluster
Note: uses the standard MatLab function “clusterdata”
4.	Assigning clusters to anatomical areas using a winner-take-all principle
MatLab: run_dwarfs_winner
5.	Stability of functional parcellations
MatLab: run_dwarfs_DCBC 
6.	Similarity of functional parcellations
MatLab: run_dwarfs_rand
Note: uses the publicly available rand MatLab function from:
https://nl.mathworks.com/matlabcentral/fileexchange/49908-adjusted-rand-index
7.	Cross-validation: split-half reliability (ALE analyses of random halves of the data)
MatLab: run_dwarfs_validate
8.	Cross-validation: split-half reliability (Pearson and Rand indices of winning clusters)
MatLab: run_dwarfs_val_winner
Note: uses the same MatLab functions described above

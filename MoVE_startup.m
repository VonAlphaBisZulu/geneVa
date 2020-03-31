%%% Modify the paths accordingly for your installation
%%% CNA: https://www2.mpi-magdeburg.mpg.de/projects/cna/download.html
%%% Ensure paths are also set correctly in CellNetAnalyzer/startcna.m lines 5 and 6
cna_root = ''; % e.g. 'CellNetAnalyzer'
geneVa_path = fullfile(pwd(), 'geneVa');

addpath(genpath(geneVa_path));

cd(cna_root);
startcna(1);

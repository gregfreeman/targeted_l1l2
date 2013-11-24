% function setup_targeted_l1l2
function setup_targeted_l1l2 

fullpath = mfilename('fullpath');
idx      = find(fullpath == filesep);
proj_folder = fullpath(1:(idx(end)-1));

addpath(proj_folder )
addpath(fullfile(proj_folder,'image_utility_toolbox'))
addpath(fullfile(proj_folder,'gf_utility_toolbox'))
addpath(fullfile(proj_folder,'experiment_framework'))
addpath(fullfile(proj_folder,'spotbox'))
addpath(fullfile(proj_folder,'compressive_sensed_image'))
addpath(fullfile(proj_folder,'image_quality_toolbox'))
addpath(fullfile(proj_folder,'libsvm'))
addpath(fullfile(proj_folder,'gf_fitting_toolbox'))
addpath(fullfile(proj_folder,'gf_solvers'))
addpath(fullfile(proj_folder,'wavelab850'))
addpath(fullfile(proj_folder,'tfocs'))
addpath(fullfile(proj_folder,'compressive_sensed_image/reconstruct'))
addpath(fullfile(proj_folder,'compressive_sensed_image/structure'))
addpath(fullfile(proj_folder,'compressive_sensed_image/reconstruct/SparseLab/SharedTools'))
addpath(fullfile(proj_folder,'compressive_sensed_image/reconstruct/SparseLab/Transforms'))
addpath(fullfile(proj_folder,'libsvm/matlab'))
addpath(fullfile(proj_folder,'wavelab850/Orthogonal'))


setup_experiment_framework


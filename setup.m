function setup()
% Run all setup functions in their respective subfolders
%
% This will create .mat files. This only has to be performed once on a
% system.
rootdir = pwd;

cd([rootdir, filesep, 'options'])
gen_param_defaults;

cd([rootdir, filesep, 'options', filesep,'interpolation'])
gen_interpolation;

cd([rootdir, filesep,'options', filesep, 'mesh'])
gen_mesh;

cd([rootdir, filesep,'options', filesep,'static_uniform'])
gen_static_uniform;

cd([rootdir, filesep,'options', filesep,'new_longer'])
gen_new_longer;

cd([rootdir, filesep, 'experiments'])
gen_def_expt_options;

cd([rootdir, filesep, 'experiments', filesep, 'alan'])
gen_alan;

cd([rootdir, filesep, 'experiments', filesep, 'barry'])
gen_barry;

cd([rootdir, filesep, 'experiments', filesep, 'carly'])
gen_carly;

cd(rootdir)
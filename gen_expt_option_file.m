expt_option_filename = 'experiments/def_expt_options.mat';

% Program to run for experiments
program_name = 'burg2';
%program_name = 'burg2';
%program_name = 'burgersTom';

% Name of parameter file to use
param_file = 'only_cubic_lagrange.mat';

% Name of the experiment output
out_filename = 'experiments/cubic_lagrange.mat';

bigNX = [20:1:320 , 325:5:640 , 640:10:1000];
%bigNT = [40,80,160]';
bigNT = [40,80,160]';

save(expt_option_filename, 'program_name', 'param_file', 'out_filename',...
    'bigNX', 'bigNT');
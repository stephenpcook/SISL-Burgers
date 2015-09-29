expt_option_filename = 'experiments/test.mat';

% Program to run for experiments
program_name = 'burg2';
%program_name = 'burg2';
%program_name = 'burgersTom';

% Name of parameter file to use
param_file = 'tmp.mat';

% Name of the experiment output
out_filename = 'experiments/test.mat';

bigNX = 20:20:320;
bigNT = [40,80,160]';

save(expt_option_filename, 'program_name', 'param_file', 'out_filename',...
    'bigNX', 'bigNT');
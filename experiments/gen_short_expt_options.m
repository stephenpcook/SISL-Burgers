bigNT = [40;60;80];
bigNX = 20:1:50;
out_filename = 'experiments/short_test.mat';
param_file = 'tmp.mat';
program_name = 'burg2';
save('short_expt_options.mat', 'bigNT', 'bigNX',...
    'out_filename','param_file','program_name');
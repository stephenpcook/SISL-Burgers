bigNT = [40;80;160];
bigNX = [(20:320),(325:5:640),(650:10:1000)];
out_filename = 'experiments/test_out.mat';
param_file = 'tmp.mat';
program_name = 'burg2';

save('def_expt_options.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');
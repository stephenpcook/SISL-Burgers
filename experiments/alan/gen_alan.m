clear all
bigNT = [40;80;160];
bigNX = [(20:320),(325:5:640),(650:10:1000)];
program_name = 'burgersSLMM';

out_filename = 'experiments/alan/alan1_out.mat';
param_file = 'alpha25.mat';
save('alan1.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

out_filename = 'experiments/alan/alan2_out.mat';
param_file = 'alpha5.mat';
save('alan2.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');


out_filename = 'experiments/alan/alan3_out.mat';
param_file = 'theta51.mat';
save('alan3.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

out_filename = 'experiments/alan/alan4_out.mat';
param_file = 'moving_exact_arc.mat';
save('alan4.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

out_filename = 'experiments/alan/alan5_out.mat';
param_file = 'only_cubic_lagrange.mat';
save('alan5.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

clear all
bigNT = [40;80;160];
bigNX = [(20:320),(325:5:640),(650:10:1000)];
program_name = 'burgersSLMM';

%% alan0
% Defaults
out_filename = 'experiments/alan/alan0_out.mat';
save('alan0.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'program_name');

%% alan1
% alpha = 0.25
out_filename = 'experiments/alan/alan1_out.mat';
param_file = 'alpha25.mat';
save('alan1.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% alan2
% alpha = 0.5
out_filename = 'experiments/alan/alan2_out.mat';
param_file = 'alpha5.mat';
save('alan2.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% alan3
% thetaX = thetaU = 0.55
out_filename = 'experiments/alan/alan3_out.mat';
param_file = 'theta55.mat';
save('alan3.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% alan4
% moving mesh, m = sqrt(0.1 + u_x^2)
out_filename = 'experiments/alan/alan4_out.mat';
param_file = 'moving_exact_arc.mat';
save('alan4.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% alan5
% Clagrange + limiter
out_filename = 'experiments/alan/alan5_out.mat';
param_file = 'only_cubic_lagrange.mat';
save('alan5.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% alan6
% Clagrange (no limiter)
out_filename = 'experiments/barry/alan6_out.mat';
param_file = 'cubic_lagrange_no_limiter.mat';

save('alan6.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% alan7
% moving mesh, m = sqrt(0.1 + u_xx^2)
out_filename = 'experiments/alan/alan7_out.mat';
param_file = 'moving_exact_curv.mat';
save('alan7.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

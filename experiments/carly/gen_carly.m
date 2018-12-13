function gen_carly()
bigNT = [50, 100, 200];
bigNX = [(10:1000)];
program_name = 'burgersSLMM';

%% carly0
% Description
out_filename = 'experiments/carly/carly0_out.mat';
param_file = 'new_longer_static.mat';
save('carly0.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% carly1
% Description
out_filename = 'experiments/carly/carly1_out.mat';
param_file = 'new_longer_static_hermite.mat';
save('carly1.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% carly2
% Description
out_filename = 'experiments/carly/carly2_out.mat';
param_file = 'new_longer_moving.mat';
save('carly2.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% carly_test
% moving mesh, m = sqrt(0.1 + u_xx^2)
bigNT = [50;100;200];
bigNX = [(10:10:500)];
program_name = 'burgersSLMM';

out_filename = 'experiments/carly/carly_test_out.mat';
param_file = 'new_longer_static_hermite.mat';
save('carly_test.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');
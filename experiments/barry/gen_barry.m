function gen_barry()
bigNT = [40;80;160];
bigNX = [20:1000];
program_name = 'burgersSLMM';


%% barry1
% Hermite, Hyman derivatives limited to monotonic
out_filename = 'experiments/barry/barry1_out.mat';
param_file = 'hermite.mat';
save('barry1.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% barry2
% Hermite, Hyman derivatives, not limited
out_filename = 'experiments/barry/barry2_out.mat';
param_file = 'hermite_no_limiter.mat';
save('barry2.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% barry3
% ENO, flux limiter
out_filename = 'experiments/barry/barry3_out.mat';
param_file = 'ENO.mat';
save('barry3.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');

%% barry4
% ENO, no flux limiter
out_filename = 'experiments/barry/barry4_out.mat';
param_file = 'ENO_no_limiter.mat';
save('barry4.mat',...
    'bigNT', 'bigNX', 'out_filename',...
    'param_file', 'program_name');


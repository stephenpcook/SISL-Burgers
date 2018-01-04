function gen_static_uniform()
% Script for generating .mat files in /static_uniform
%
% Creates alpha25.mat and alpha5.mat
x_l = -1;
x_r = 4;
c = 1;
t0 = 0;
epsilon = 0.0001;

%% alpha25

alpha_0 = 0.25;

u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
u_l = u0(x_l);
u_r = u0(x_r);
save('alpha25','alpha_0','u0','u_l','u_r');

%% alpha5
alpha_0 = 0.5;

u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
u_l = u0(x_l);
u_r = u0(x_r);
save('alpha5','alpha_0','u0','u_l','u_r');

%% theta55
theta_t = 0.55;
theta_x = 0.55;
save('theta55','theta_t','theta_x');

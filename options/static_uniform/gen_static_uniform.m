% Script for generating .mat files in /static_uniform
%
% Creates alpha25.mat and alpha5.mat
clear all
alpha_0 = 0.25;
save('alpha25','alpha_0');

clear all
alpha_0 = 0.5;
save('alpha5','alpha_0');

clear all
theta_t = 0.51;
theta_x = 0.51;
save('theta51','theta_t','theta_x');

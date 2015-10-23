% Script for generating .mat files in /static_uniform
%
% Creates alpha25.mat and alpha5.mat
clear all
alpha_0 = 0.25;
save('alpha25','alpha_0');

clear all
alpha_0 = 0.5;
save('alpha5','alpha_0');
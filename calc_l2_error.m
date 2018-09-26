function l2_error = calc_l2_error(X, U, param_file)
%CALC_L2_ERROR Summary of this function goes here
%   Detailed explanation goes here
old_path = addpath([pwd,'\interpolation'],...
                   [pwd,'\mm_suite'],...
                   [pwd,'\options']);

load('params_default');
load(param_file);

X_full = [x_l; X; x_r];
U_full = [u_l; U; u_r];

u_final = @(x) u0(x - c*tmax);

l2_error = sqrt(trapz(X_full, (U_full - u_final(X_full)).^2));

end


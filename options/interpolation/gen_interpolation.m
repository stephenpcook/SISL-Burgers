clear all
%% only_cubic_lagrange
interpolation = 'CLagrange';
limiter = 1;
save('only_cubic_lagrange',...
    'interpolation', 'limiter');

%% cubic_lagrange_no_limiter
interpolation = 'CLagrange';
limiter = 0;
save('cubic_lagrange_no_limiter',...
    'interpolation', 'limiter');
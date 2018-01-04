function gen_interpolation()
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

%% hermite
interpolation = 'hermite';
limiter = 1;
save('hermite',...
    'interpolation','limiter');


%% hermite_no_limiter
interpolation = 'hermite';
limiter = 0;
save('hermite_no_limiter',...
    'interpolation','limiter');

%% ENO
interpolation = 'ENO';
limiter = 1;
save('ENO',...
    'interpolation','limiter');

%% ENO_no_limiter
interpolation = 'ENO';
limiter = 0;
save('ENO_no_limiter',...
    'interpolation','limiter');

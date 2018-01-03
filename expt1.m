function out_filename = expt1(expt_option_file)
%EXPT1 Runs a given SISL code a number of times%
%
% out_filename = EXPT1(expt_option_file) runs one of the programs for a
% number of NX and NT, as defined in the experiment file expt_option_file
% and saves the results, then runs expt2.
%
% See also: EXPT2 GEN_PARAM_DEFAULTS

%  IN:
%   expt_option_file - string with filename.mat
%
%  option_file should define
%    program_name
%    out_filename
%    bigNX
%    bigNT
%    param_file
%
%  param_file is either the name of the burgers option file to use, or just
% use the same one for everything (so put the above things in param_file).

%%%%%%%%%%%%%%
% Load files %
%%%%%%%%%%%%%%
% % These to be loaded from expt_option_file
%program_name = 'burg2';
%out_filename = 'test.mat';
%bigNX = 20:20:320;
%bigNT = [40,80,160]';
% % Optional: If not present use expt_option_file
%param_file = 'tmp.mat';

addpath('.');
addpath(genpath('experiments'));
addpath(genpath('options'));

load('params_default','c','alpha_0','epsilon','tmax')
load(expt_option_file , 'program_name', 'param_file', 'out_filename',...
    'bigNX', 'bigNT')

if isempty(who('param_file'))
    % If param_file is not defined in expt_option_file, then pass
    % expt_option_file to the code for options.
    param_file = expt_option_file;
else
    load(param_file); %#ok<NODEF> Warning check in if statement
end % if

if isempty(who('out_filename'))
    % Default for out_filename
    out_filename = 'tmp.mat';
end % if
%%%%%%%%%%%%%%%%%%
% Initialisation %
%%%%%%%%%%%%%%%%%%
% grad is the gradient of the exact solution of a tanh profile.
%grad = -alpha_0^2/(2*epsilon);

bigC = zeros(length(bigNX),length(bigNT));
bigC2 = zeros(length(bigNX),length(bigNT));
bigM = zeros(length(bigNX),length(bigNT));
bigW = zeros(length(bigNX),length(bigNT));
bigMinDx = zeros(length(bigNX),length(bigNT));

%%%%%%%%%%%%%
% Main code %
%%%%%%%%%%%%%
for ii = 1:length(bigNX)
  for jj = 1:length(bigNT)
    N = bigNX(ii);
    tN = bigNT(jj);
    switch program_name
      case 'burgersSLMM'
        %[U,X] = burgersSLMM(N,tN,param_file);
        [U,X,X_star,DxMin] = burgersSLMM(N,tN,param_file);
      case 'burg2'
        [U,X,X_star] = burg2(N,tN,param_file);
        DxMin = min(diff(X));
      case 'burgersTom'
        [U,X] = burgersTom(N,tN);
        DxMin = min(diff(X));
      otherwise
        exec(['[U,X] = ',program_name,'(N,tN);']);
    end
    [m,x_star] = get_m_x(U,X,c,alpha_0);
    bigC(ii,jj) = x_star/tmax;
    % Alternative way of calculating speed: line of best fit
    T = linspace(0,1.5,tN+1)';
    f = fit(T,X_star,'poly1');
    bigC2(ii,jj) = f.p1;
    %
    bigM(ii,jj) = m;
    bigMinDx(ii,jj) = min(DxMin);
  end % for jj
end % for ii

[BigNT, BigNX] = meshgrid(bigNT,bigNX);  %#ok<ASGLU> Warning unused; save

% bigEps extracted from the midpoint gradient
bigEps = -0.5*alpha_0^2./bigM;  %#ok<*NASGU> File warning unused; save

%display(bigNX)
%display(bigNT)
%display(bigC)

%display(grad)
%display(bigEps)

%%%%%%%%%%%%%%%%%%%%%
% Save to .mat file %
%%%%%%%%%%%%%%%%%%%%%
save(out_filename, 'param_file', 'program_name',...
    'c','epsilon','alpha_0','tmax',...
    'bigC','bigC2','bigEps','bigM','bigMinDx',...
  'bigNX','bigNT','BigNX','BigNT');

fprintf(['Values printed to ',out_filename,'\n'])

expt2(out_filename)

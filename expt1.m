function out_filename = expt1(expt_option_file)
% Runs a given SISL code a number of times, defined by expt_option_file
%
% Runs one of the programs for a number of NX and NT, and saves and
% tabulates the results.
%
% function out_filename = expt1(expt_option_file)

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
bigM = zeros(length(bigNX),length(bigNT));
bigW = zeros(length(bigNX),length(bigNT));

%%%%%%%%%%%%%
% Main code %
%%%%%%%%%%%%%
for ii = 1:length(bigNX)
  for jj = 1:length(bigNT)
    N = bigNX(ii);
    tN = bigNT(jj);
    switch program_name
      case 'burgersSLMM'
        [U,X] = burgersSLMM(N,tN,param_file);
      case 'burg2'
        [U,X] = burg2(N,tN,param_file);
      case 'burgersTom'
        [U,X] = burgersTom(N,tN);
      otherwise
        exec(['[U,X] = ',program_name,'(N,tN);']);
    end
    [m,x_star,w] = get_m_x(U,X,c,alpha_0);
    bigC(ii,jj) = x_star/tmax;
    bigM(ii,jj) = m;
    bigW(ii,jj) = w;
  end % for jj
end % for ii

[BigNT, BigNX] = meshgrid(bigNT,bigNX);  %#ok<ASGLU> Warning unused; save

% bigEps1 extracted from the midpoint gradient
bigEps1 = -0.5*alpha_0^2./bigM;  %#ok<*NASGU> File warning unused; save
% bigEps2 extracted from the front width
bigEps2 = alpha_0*bigW/(4*atanh(0.95)); 

%display(bigNX)
%display(bigNT)
%display(bigC)

%display(grad)
%display(bigEps1)
%display(bigEps2)

%%%%%%%%%%%%%%%%%%%%%%%
% Create LaTeX Tables %
%%%%%%%%%%%%%%%%%%%%%%%
mystr = mk_latex_table(bigC,bigNX,bigNT,'Estimate of C','%.3f',1,0);
mystr2 = mk_latex_table(-0.5*alpha_0^2./bigM,bigNX,bigNT,...
    'Estimate of $\\eps$ from gradient at $x=c$.','%.5f',1,0);
mystr3 = mk_latex_table(bigW.*(alpha_0/(4*1.832)),bigNX,bigNT,...
    'Estimate of $\\eps$ from width of front.','%.5f',1,0);

%%%%%%%%%%%%%%%%%%%%%
% Save to .mat file %
%%%%%%%%%%%%%%%%%%%%%
save(out_filename, 'c','epsilon','alpha_0','tmax',...
    'bigC','bigEps1','bigEps2','bigM',...
  'bigNX','bigNT','BigNX','BigNT','mystr','mystr2');

fprintf(['Values printed to ',out_filename,'\n'])

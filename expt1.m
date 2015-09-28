% Runs one of the programs for a number of N and TN, and tabulates the
% results.

% TODO - make this a function which takes
%  in:
%   expt_option_file - string with filename.mat
%
%  option_file should define 
%    program_name
%    out_filename
%    bigNX
%    bigNT
%    ~param_file~
%
%  param_file is either the name of the burgers option file to use, or just 
% use the same one for everything (so put the above things in param_file).

load('params_default','c','alpha_0','epsilon','tmax')
%load(expt_option_file)

% Get rid of this
plotting = 0;
save('tmp.mat','plotting')
param_file = 'tmp.mat';
%

% These to be loaded, get rid of this
program_name = 'burg2';
%program_name = 'burg2';
%program_name = 'burgersTom';
out_filename = 'test.mat';
bigNX = 20:20:320;
bigNT = [40,80,160]';
%

grad = -alpha_0^2/(2*epsilon);

bigC = zeros(length(bigNX),length(bigNT));
bigM = zeros(length(bigNX),length(bigNT));
bigM2 = zeros(length(bigNX),length(bigNT));


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
    %bigM2(ii,jj) = -m2;
  end % for jj
end % for ii

[BigNT, BigNX] = meshgrid(bigNT,bigNX);

bigEps1 = -0.5*alpha_0^2./bigM;
bigEps2 = alpha_0*bigW/(4*atanh(0.95));

%display(bigNX)
%display(bigNT)
%display(bigC)

%display(grad)
%display(bigEps1)
%display(bigEps2)

mystr = mk_latex_table(bigC,bigNX,bigNT,'Estimate of C','%.3f',1,0);
mystr2 = mk_latex_table(-0.5*alpha_0^2./bigM,bigNX,bigNT,...
    'Estimate of $\\eps$ from gradient at $x=c$.','%.5f',1,0);
mystr3 = mk_latex_table(bigW.*(alpha_0/(4*1.832)),bigNX,bigNT,...
    'Estimate of $\\eps$ from width of front.','%.5f',1,0);

save(out_filename, 'c','epsilon','alpha_0','tmax',...
    'bigC','bigEps1','bigEps2','bigM',...
  'bigNX','bigNT','BigNX','BigNT','mystr','mystr2');


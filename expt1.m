% Runs one of the programs for a number of N and TN, and tabulates the
% results.

program_name = 'burg2';
%program_name = 'burg2';
%program_name = 'burgersTom';
filename_end = '_high_c.mat';
filename = ['out_',program_name,filename_end];
%bigN = [20,40,80,160,320,640];
%bigTN = [20,40,80,160,320,640]';
%bigN = [20,40,80,160,320,640];
%bigTN = [20,40,80,160,320,640]';
bigN = 20:1:640;
bigTN = [40,80,160]';

c = 1;
alpha_0 = 0.1;
epsilon=0.0001;
t_max = 1.5;
grad = -alpha_0^2/(2*epsilon);

bigC = zeros(length(bigN),length(bigTN));
bigM = zeros(length(bigN),length(bigTN));
bigM2 = zeros(length(bigN),length(bigTN));


for ii = 1:length(bigN)
  for jj = 1:length(bigTN)
    N = bigN(ii);
    tN = bigTN(jj);
    switch program_name
      case 'burgersSLMM'
        [U,X] = burgersSLMM(N,tN);
      case 'burg2'
        [U,X] = burg2(N,tN);
      case 'burgersTom'
        [U,X] = burgersTom(N,tN);
      otherwise
        exec(['[U,X] = ',program_name,'(N,tN);']);
    end
    [m,cDt,w] = get_m_x(U,X,c,alpha_0);
    bigC(ii,jj) = cDt/t_max;
    bigM(ii,jj) = m;
    bigW(ii,jj) = w;
    %bigM2(ii,jj) = -m2;
  end % for jj
end % for ii

[BigTN, BigN] = meshgrid(bigTN,bigN);

bigEps1 = -0.5*alpha_0^2./bigM;
bigEps2 = alpha_0*bigW/(4*atanh(0.95));

display(bigN)
display(bigTN)
display(bigC)

display(grad)
display(bigEps1)
display(bigEps2)
mystr = mk_latex_table(bigC,bigN,bigTN,'Estimate of C','%.3f',1,0);
mystr2 = mk_latex_table(-0.5*alpha_0^2./bigM,bigN,bigTN,...
    'Estimate of $\\eps$ from gradient at $x=c$.','%.5f',1,0);
mystr3 = mk_latex_table(bigW.*(alpha_0/(4*1.832)),bigN,bigTN,...
    'Estimate of $\\eps$ from width of front.','%.5f',1,0);
save(filename, 'c','epsilon','t_max','bigC','bigEps1','bigEps2','bigM',...
  'bigN','bigTN','BigN','BigTN','mystr','mystr2');


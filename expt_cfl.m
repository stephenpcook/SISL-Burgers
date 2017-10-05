function expt_cfl
% EXPT_CFL
%
%
load('options\params_default.mat')
load('options\interpolation\hermite.mat')
c = 1;
alpha_0 = 0.1;
tmax = 1.5;

CFL_vec = [0.7, 1.4, 2.1];
%CFL_vec = [2.1];
Nx_vec = 10:1000;
L = x_r - x_l;
%tmax;
%c;

% Inputs
[bigNx,CFL] = meshgrid(Nx_vec,CFL_vec);
% Outputs
bigEps = zeros(size(bigNx));
bigC = bigEps;

bigNt = floor((tmax .* c .* bigNx) ./(CFL * L));
for ii = 1:length(CFL_vec)
  tic
  fprintf('CFL number %f ... ',CFL_vec(ii));
  parfor jj=1:length(Nx_vec)
    Nx = bigNx(ii,jj);
    Nt = bigNt(ii,jj);
    [U,X,X_star] = burgersSLMM(Nx,Nt,'options\params_default.mat');
    % Calculate the viscosity parameter at end time
    [m,~] = get_m_x(U,X,c,alpha_0);
    bigEps(ii,jj) = -0.5*alpha_0^2./m;
    % Calculate the average wave speed
    T = linspace(0,tmax,Nt+1)';
    if (Nt==1)
        X_star = X_star.';
    end
    f = fit(T,X_star,'poly1');
    bigC(ii,jj) = f.p1;
  end % for jj
  fprintf('done.\n');
  toc
  save('CFL_WIP.txt','CFL_vec','bigNx','bigC','bigEps','-ascii')
  save('CFL_WIP.mat','CFL_vec','bigNx','bigC','bigEps')
end % for ii

save('CFL4.txt','CFL_vec','bigNx','bigC','bigEps','-ascii')
save('CFL4.mat','CFL_vec','bigNx','bigC','bigEps')

expt_cfl_plot('CFL4.mat');
end % function expt_cfl

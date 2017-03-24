function [h1] = theory_plots(in_filename)
% theory_plots
if nargin==0
    in_filename = 'test.mat';
end

load(in_filename)

if ~exist('BigNX','var')
    [BigNT, BigNX] = meshgrid(bigNT,bigNX);
end
if ~exist('c','var')
    c = 1;
end
if ~exist('alpha_0','var')
    alpha_0 = 0.1;
end
if ~exist('epsilon','var')
    epsilon = 1e-4;
end

Dx = 5./BigNX;
Dt = 1.5./BigNT;

BigEll = floor(c*Dt./Dx);
epsmin = alpha_0*Dx./4;

% Equation eq:eps_hat in PaperNew.tex
BigEpsHat = epsilon + 0.5*((2*BigEll+1).*Dx - c*Dt)*c...
    - (BigEll.^2 + BigEll).*(Dx.^2./(2*Dt));
BigEpsTheo = max(BigEpsHat, epsmin);

% Equation eq:c_hat_theory in PaperNew.tex
BigCHatTheo = c - alpha_0^2./(6*BigEpsTheo)...
    .*((2*BigEll+1).*Dx - 2*c*Dt);



% Plotting
set(groot,'DefaultTextInterpreter','latex');
%set(groot,'DefaultLegendInterpreter','latex');
% Set linecolor to default black
set(groot,'defaultAxesColorOrder',[0 0 0]);
% Set lines to cycle these styles
myLines = {'-','--','-.'};
set(groot,'defaultAxesLineStyleOrder',myLines)

for ii = 1:length(bigNT)
  myLegend{ii} = ['N_t = ',num2str(bigNT(ii))];
end % for ii



h1 = figure;
semilogy(BigNX, BigEpsTheo);
xlabel('$N_x$')
ylabel('$\widehat{\varepsilon}$')
legend(myLegend)



h2 = figure;
plot(BigNX, BigCHatTheo);
xlabel('$N_x$')
ylabel('$\widehat{c}$')
ylim([0.9,1.1])
legend(myLegend)

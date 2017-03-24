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

h1 = figure;
semilogy(BigNX, BigEpsTheo);
xlabel('N_x')
ylabel('\epsilon_{hat}')
%title('Theoretical \epsilon_{hat}')

h2 = figure;
plot(BigNX, BigCHatTheo);
xlabel('N_x')
ylabel('\c_hat')
ylim([0.9,1.1])
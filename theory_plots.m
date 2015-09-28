function [h1] = theory_plots(in_filename)
% theory_plots
if nargin==0
    in_filename = 'test.mat';
end

load(in_filename)

Dx = 5./BigNX;
Dt = 1.5./BigNT;

BigEll = floor(c*Dt./Dx);
epsmin = alpha_0*Dx./4;

BigEpsHat = epsilon + 0.5*((2*BigEll+1).*Dx - Dt)...
    - (BigEll.^2 + BigEll).*(Dx.^2./(2*Dt));
BigEpsTheo = max(BigEpsHat, epsmin);

h1 = figure;
loglog(BigNX, BigEpsTheo);
xlabel('N_x')
ylabel('\epsilon_{hat}')
title('Theoretical \epsilon_{hat}')
function [] = expt2(in_filename)
%EXPT2 Plots experiment results created by expt1
%
% EXPT2(in_filename) creates black and white plots of epsilon and k created
% by expt1. Can be modified to create plots of the minimum mesh spacing.
%
% See also: EXPT1

%if ( nargin==0 )
%  in_filename = 'experiments/test_out.mat';
%end % if nargin
load(in_filename)


for ii = 1:length(bigNT)
  myLegend{ii} = ['N_t = ',num2str(bigNT(ii))];
end % for ii

% Set linecolor to default black
set(groot,'defaultAxesColorOrder',[0 0 0]);
% Set lines to cycle these styles
myLines = {'-','--','-.'};
set(groot,'defaultAxesLineStyleOrder',myLines);



h1 = figure;
set(h1,'defaulttextinterpreter','latex');
loglog(bigNX,bigEps(:,1),myLines{1})
hold on
for ii=2:length(bigNT)-1
loglog(bigNX,bigEps(:,ii),myLines{2})
end % for ii
loglog(bigNX,bigEps(:,end),myLines{3})
title([program_name,' : bigEps for different $N_t$'])
xlabel('$N_x$')
ylabel('$\widehat{\varepsilon}$')
legend(myLegend)


h3 = figure;
set(h3,'defaulttextinterpreter','latex');
plot(bigNX, bigC2)
xlabel '$N_x$'
ylabel '$k$'
legend(myLegend)
%title('Estimate of k from average speed of front (grad of line of best fit for position.')

% figure
% loglog(bigNX, bigMinDx)
% xlabel 'N_x'
% ylabel 'min(\Delta X)'
% legend(myLegend)

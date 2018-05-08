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
set(groot,'defaultAxesColorOrder',...
    [0 0 0;
     0.4 0.4 0.4
     0.6 0.6 0.6]);
% Set lines to cycle these styles
myLines = {'-','--','-.'};
set(groot,'defaultAxesLineStyleOrder',myLines);



h1 = figure;
set(h1,'defaulttextinterpreter','latex');
semilogy(bigNX,bigEps(:,1),myLines{1})
hold on
for ii=2:length(bigNT)-1
semilogy(bigNX,bigEps(:,ii),myLines{2})
end % for ii
semilogy(bigNX,bigEps(:,end),myLines{3})
%title([program_name,' : bigEps for different $N_t$'])
xlabel('$N_x$')
ylabel('$\widehat{\varepsilon}$')
legend(myLegend)
ylim([1e-4, 1e-1])


h3 = figure;
set(h3,'defaulttextinterpreter','latex');
plot(bigNX, bigC2(:,1))
hold on
for ii=2:length(bigNT)-1
plot(bigNX,bigC2(:,ii),myLines{2})
end % for ii
semilogy(bigNX,bigC2(:,end),myLines{3})
hold off
xlabel '$N_x$'
ylabel '$\widehat{c}$'
legend(myLegend)
ylim([c-alpha_0, c+alpha_0])
%title('Estimate of k from average speed of front (grad of line of best fit for position.')

h4 = figure;
set(h4,'defaulttextinterpreter','latex');
semilogy(bigNX, bigMinDx(:,1),myLines{1})
hold on
for ii=2:length(bigNT)-1
semilogy(bigNX, bigMinDx(:,ii),myLines{2})
end % for ii
semilogy(bigNX, bigMinDx(:,end),myLines{3})
hold off
xlabel('$N_x$','FontSize',18)
ylabel('$\displaystyle{\min_{j,n}(X^n_j - X^n_{j-1})}$','FontSize',18)
legend(myLegend,'FontSize',18)
ylim([1e-5-eps 1])

function [] = expt2(in_filename)
%EXPT2 Plots experiment results created by expt1
%
% EXPT2(in_filename) creates black and white plots of epsilon and k created
% by expt1. Can be modified to create plots of the minimum mesh spacing.
%
% See also: EXPT1


load(in_filename)

for ii = 1:length(bigNT)
  myLegend{ii} = ['\Delta t = ',num2str(tmax / bigNT(ii))];
end % for ii

% Set lines to cycle these styles and colors
myLines = {'-','--','-.'};
myColors = {[0 0 0], [0.4 0.4 0.4], [0.6 0.6 0.6]};
FontSize = 18;
LegendFontSize = 18;

h1 = figure;
set(h1,'defaulttextinterpreter','latex');
a = semilogy(bigNX, bigEps(:,1));
set(a, 'Color', myColors{1}, 'LineStyle', myLines{1});
hold on
for ii=2:length(bigNT)-1
  a = semilogy(bigNX, bigEps(:,ii));
set(a, 'Color', myColors{2}, 'LineStyle', myLines{2});
end % for ii
a = semilogy(bigNX,bigEps(:,end));
set(a, 'Color', myColors{end}, 'LineStyle', myLines{end});
%title([program_name,' : bigEps for different $N_t$'])
xlabel('$N_x$', 'FontSize', FontSize)
ylabel('$\widehat{\varepsilon}$', 'FontSize', FontSize)
legend(myLegend)
ylim([1e-4, 1e-1])


h3 = figure;
set(h3,'defaulttextinterpreter','latex');
a = plot(bigNX, bigC2(:,1));
set(a, 'Color', myColors{1}, 'LineStyle', myLines{1});
hold on
for ii=2:length(bigNT)-1
  a = plot(bigNX, bigC2(:,ii));
  set(a, 'Color', myColors{2}, 'LineStyle', myLines{2});
end % for ii
a = semilogy(bigNX, bigC2(:,end));
set(a, 'Color', myColors{3}, 'LineStyle', myLines{3});
hold off
xlabel('$N_x$', 'FontSize', FontSize)
ylabel('$\widehat{c}$', 'FontSize', FontSize)
legend(myLegend)
ylim([c-alpha_0, c+alpha_0])
%title('Estimate of k from average speed of front (grad of line of best fit for position.')

h4 = figure;
set(h4,'defaulttextinterpreter','latex');
a = semilogy(bigNX, bigMinDx(:,1));
set(a, 'Color', myColors{1}, 'LineStyle', myLines{1});
hold on
for ii=2:length(bigNT)-1
  a = semilogy(bigNX, bigMinDx(:,ii));
  set(a, 'Color', myColors{2}, 'LineStyle', myLines{2});
end % for ii
a = semilogy(bigNX, bigMinDx(:,end));
set(a, 'Color', myColors{3}, 'LineStyle', myLines{3});
hold off
xlabel('$N_x$', 'FontSize', FontSize)
ylabel('$\displaystyle{\min_{j,n}(X^n_j - X^n_{j-1})}$', 'FontSize', FontSize)
legend(myLegend, 'FontSize', LegendFontSize)
ylim([1e-5-eps 1])

h5 = figure;
set(h5,'defaulttextinterpreter','latex');
set(h5,'defaultAxesColorOrder',...
    [0 0 0;
     0.4 0.4 0.4
     0.6 0.6 0.6]);
a = semilogy(bigNX, bigL2(:,1));
set(a, 'Color', myColors{1}, 'LineStyle', myLines{1});
hold on
for ii=2:length(bigNT)-1
  a = semilogy(bigNX, bigL2(:,ii));
  set(a, 'Color', myColors{2}, 'LineStyle', myLines{2});
end % for ii
a = semilogy(bigNX, bigL2(:,end));
set(a, 'Color', myColors{3}, 'LineStyle', myLines{3});
hold off
xlabel('$N_x$', 'FontSize', FontSize)
ylabel('$\displaystyle{||U - u(x, 1.5)||_2}$', 'FontSize', FontSize)
ylim([1e-3, 1e-0])
legend(myLegend, 'FontSize', LegendFontSize)

function expt_cfl_plot(cfl_mat_file)
if nargin<1
    cfl_mat_file = 'CFL.mat';
end % if nargin
load(cfl_mat_file,'CFL_vec','bigNx','bigC','bigEps')

% Set linecolor to default black
set(groot,'defaultAxesColorOrder',[0 0 0]);
% Set lines to cycle these styles
myLines = {'-','--','-.'};
set(groot,'defaultAxesLineStyleOrder',myLines);

for ii = 1:length(CFL_vec)
    myLegend{ii} = ['CFL = ',num2str(CFL_vec(ii))];
end

h1 = figure(1);
set(h1,'defaulttextinterpreter','latex')
plot(bigNx',bigC')
ylim([0.9,1.1])
xlabel('$N_x$','FontSize',18)
ylabel('$\widehat{c}$','FontSize',18)
legend(myLegend,'FontSize',18)

h2 = figure(2);
set(h2,'defaulttextinterpreter','latex')
semilogy(bigNx',bigEps')
xlabel('$N_x$','FontSize',18)
ylabel('$\widehat{\varepsilon}$','FontSize',18)
legend(myLegend,'FontSize',18)

end
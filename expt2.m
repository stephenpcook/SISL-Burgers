function [] = expt2(program_name,in_filename)
% TODO - Make this a function that takes the in_filename as an input.

%close all
%progs = {'burgersSLMM','burg2','burgersTom'};
%progs = {'burgersSLMM_MM'};
%program_name = 'burg2';
%in_filename = 'experiments/test_out.mat';
load(in_filename)


for ii = 1:length(bigNT)
  myLegend{ii} = ['N_t = ',num2str(bigNT(ii))];
end % for ii


% figure
% loglog(bigNX,(2/c*bigEps(:,1)+1.5*c/bigNT(1)).*bigNX','k-')
% hold on
% xlabel 'N_x'
% ylabel 'K'
% for ii = 2:length(bigNT)-1
% loglog(bigNX,(2/c*bigEps(:,ii)+1.5*c/bigNT(ii)).*bigNX','b-')
% end % for ii
% loglog(bigNX,(2/c*bigEps(:,end)+1.5*c/bigNT(end)).*bigNX','r-')
% hold off
% title([program_name,' : N*(bigEps/0.5v + vDt) for different tN'])
% xlabel 'N_x'
% ylabel 'K'
% legend(myLegend)


figure
loglog(bigNX,bigEps(:,1),'k-')
hold on
for ii=2:length(bigNT)-1
loglog(bigNX,bigEps(:,ii),'b-')
end % for ii
loglog(bigNX,bigEps(:,end),'r-')
title([program_name,' : bigEps for different tN'])
xlabel 'N_x'
ylabel 'bigEps'
legend(myLegend)

figure
plot(bigNX, bigC)
xlabel 'N_x'
ylabel 'bigC'
legend(myLegend)
title('Estimate of k from endpoint of front.')

figure
plot(bigNX, bigC2)
xlabel 'N_x'
ylabel 'bigC2'
legend(myLegend)
title('Estimate of k from average speed of front (grad of line of best fit for position.')

figure
loglog(bigNX, bigMinDx)
xlabel 'N_x'
ylabel 'min(\Delta X)'
legend(myLegend)

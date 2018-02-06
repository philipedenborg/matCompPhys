close all
data = importdata('cutoff_file');
plot(data(:,1),data(:,2),'-o')

xlabel('Cut-off energy [eV]','interpreter','latex','fontsize',16)
ylabel('Potential energy [eV]','interpreter','latex','fontsize',16)

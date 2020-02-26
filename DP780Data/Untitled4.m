data = importdata('epsc3.out');
data = data.data;

figure;
hold on;
plot(data(:,1),data(:,7));
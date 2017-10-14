clear all;
load 'exact.m';

figure(1)
plot(exact(:,1),exact(:,2));
xlabel('r');
ylabel('\rho');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);

figure(2)
plot(exact(:,1),exact(:,3));
xlabel('r');
ylabel('u');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);

figure(3)
plot(exact(:,1),exact(:,4))
xlabel('r');
ylabel('p');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);

figure(4)
plot(exact(:,1),exact(:,5))
xlabel('r');
ylabel('e');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
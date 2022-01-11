% Plot tilt and start of ump


n410 = 4;
n660 = 1;
time = 199;


load(['tilt_2021-4cm-' num2str(n410) '-' num2str(n660) '.mat'])


x = [time time];
y = [45 90];

tiltsm = smoothdata(tilt(:,4),'movmean',20);


plot(tilt(:,2),tilt(:,4),'k','LineWidth',0.5)
title(['Dislocation: 2cm ' num2str(n410) '-' num2str(n660)])
hold on
plot(x,y,'--','LineWidth',1.5)
plot(tilt(:,2),tiltsm,'LineWidth',2)
ylim([45 90]);
xlim([100 400]);
xlabel('Distance from center (km)')
ylabel('Radius (km)')
legend('Tilt','UMP start','Smoothed 5 Myr')
xlabel('Time (Myr)')
ylabel('Tilt (degree)')

%saveas(gca,['tilt_2021-dis-2cm-' num2str(n410) '-' num2str(n660) '.png']);

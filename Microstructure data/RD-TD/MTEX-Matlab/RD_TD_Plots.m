figure(1)
plot(ebsd1('aluminium'),ebsd1('aluminium').orientations,'micronbar')
hold on
plot(selected_grains1.boundary,'lineWidth',0.5)
hold on
[omega,a,b] = selected_grains1.fitEllipse;
plotEllipse(selected_grains1.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd2('aluminium'),ebsd2('aluminium').orientations,'micronbar')
hold on
plot(selected_grains2.boundary,'lineWidth',0.5)
hold on
[omega,a,b] = selected_grains2.fitEllipse;
plotEllipse(selected_grains2.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd3('aluminium'),ebsd3('aluminium').orientations,'micronbar')
hold on
plot(selected_grains3.boundary,'lineWidth',0.5)
hold on
[omega,a,b] = selected_grains3.fitEllipse;
plotEllipse(selected_grains3.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd4('aluminium'),ebsd4('aluminium').orientations,'micronbar')
hold on
plot(selected_grains4.boundary,'lineWidth',0.5)
hold on
[omega,a,b] = selected_grains4.fitEllipse;
plotEllipse(selected_grains4.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd5('aluminium'),ebsd5('aluminium').orientations,'micronbar')
hold on
plot(selected_grains5.boundary,'lineWidth',0.5)
hold on
[omega,a,b] = selected_grains5.fitEllipse;
plotEllipse(selected_grains5.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold off

ipfKey = ipfColorKey(ebsd1('aluminium'));
figure(2)
plot(ipfKey)



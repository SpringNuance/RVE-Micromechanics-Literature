figure(1)
plot(ebsd1('aluminium'),ebsd1('aluminium').orientations,'micronbar')
hold on
plot(selected_grains1.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains1.fitEllipse;
plotEllipse(selected_grains1.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd2('aluminium'),ebsd2('aluminium').orientations,'micronbar')
hold on
plot(selected_grains2.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains2.fitEllipse;
plotEllipse(selected_grains2.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd3('aluminium'),ebsd3('aluminium').orientations,'micronbar')
hold on
plot(selected_grains3.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains3.fitEllipse;
plotEllipse(selected_grains3.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd4('aluminium'),ebsd4('aluminium').orientations,'micronbar')
hold on
plot(selected_grains4.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains4.fitEllipse;
plotEllipse(selected_grains4.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd5('aluminium'),ebsd5('aluminium').orientations,'micronbar')
hold on
plot(selected_grains5.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains5.fitEllipse;
plotEllipse(selected_grains5.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd6('aluminium'),ebsd6('aluminium').orientations,'micronbar')
hold on
plot(selected_grains6.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains6.fitEllipse;
plotEllipse(selected_grains6.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd7('aluminium'),ebsd7('aluminium').orientations,'micronbar')
hold on
plot(selected_grains7.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains7.fitEllipse;
plotEllipse(selected_grains7.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd8('aluminium'),ebsd8('aluminium').orientations,'micronbar')
hold on
plot(selected_grains8.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains8.fitEllipse;
plotEllipse(selected_grains8.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd9('aluminium'),ebsd9('aluminium').orientations,'micronbar')
hold on
plot(selected_grains9.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains9.fitEllipse;
plotEllipse(selected_grains9.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd10('aluminium'),ebsd10('aluminium').orientations,'micronbar')
hold on
plot(selected_grains10.boundary,'lineWidth',0.5)
[omega,a,b] = selected_grains10.fitEllipse;
plotEllipse(selected_grains10.centroid,a,b,omega,'lineColor','w','linewidth',0.2)

hold off

ipfKey = ipfColorKey(ebsd1('aluminium'));
figure(2)
plot(ipfKey)



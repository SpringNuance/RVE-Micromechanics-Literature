figure(1)
plot(ebsd1('aluminium'),ebsd1('aluminium').orientations,'micronbar')
hold on
plot(selected_grains1.boundary,'lineWidth',0.5)
hold on
plotEllipse(selected_grains1.centroid,a1,b1,omega1,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd2('aluminium'),ebsd2('aluminium').orientations,'micronbar')
hold on
plot(selected_grains2.boundary,'lineWidth',0.5)
hold on
plotEllipse(selected_grains2.centroid,a2,b2,omega2,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd3('aluminium'),ebsd3('aluminium').orientations,'micronbar')
hold on
plot(selected_grains3.boundary,'lineWidth',0.5)
hold on
plotEllipse(selected_grains3.centroid,a3,b3,omega3,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd4('aluminium'),ebsd4('aluminium').orientations,'micronbar')
hold on
plot(selected_grains4.boundary,'lineWidth',0.5)
hold on
plotEllipse(selected_grains4.centroid,a4,b4,omega4,'lineColor','w','linewidth',0.2)

hold on
plot(ebsd5('aluminium'),ebsd5('aluminium').orientations,'micronbar')
hold on
plot(selected_grains5.boundary,'lineWidth',0.5)
hold on
plotEllipse(selected_grains5.centroid,a5,b5,omega5,'lineColor','w','linewidth',0.2)

hold off

ipfKey = ipfColorKey(ebsd1('aluminium'));
figure(2)
plot(ipfKey)

% set pole figure annotation
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'RD','TD'},...
  'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
figure(3)
plotPDF(odf,[Miller(1,0,0,cs),Miller(0,1,1,cs),Miller(1,1,1,cs)],'antipodal')
figure(4)
plot(odf,'phi2',[0,15,30,45,60,75,90]*degree,'antipodal')
figure(5)
odf.SS = specimenSymmetry('mmm');
plot(odf,'phi2',[0,15,30,45,60,75,90]*degree,'antipodal')


% Misorientation distribution function
figure(6)
plotAngleDistribution(mdf)
legend('Random grains')
title('RD-TD')
ylabel('Frequency PDF,-')
%figure(7)
%plot(omega_deg,frequency_percent,'b','LineWidth',2)
%calcMDF(grains1.meanOrientation)



%% Script for EBSD Data analysis
% Purpose: EBSD analysis with Matlab-MTEX toolbox
% Functions: Import EBSD data
%            EBSD plotting & initial analyses
%            Grain reconstruction
%            Grain Size & Shape analysis
%            Grain Orientation analysis
%            Grain Misorientation analysis
%            GND analysis
%
% Record of revision:
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  26-3-2023        Rongfei Juan       Original code for MTEX 5.8.1
%
%% Import Script for EBSD Data
% Specify Crystal and Specimen Symmetries
% crystal symmetry
CS = {... 
  'notIndexed',...
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.6 0.8 1])};
% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
% path to files
pname = 'C:\Users\juanr1\OneDrive - Aalto University\Journal CP parameter calibration\Dislocation density\MTEX analysis\01 Exp data';
% which files to be imported
fname = [pname '\B_1.ctf'];
%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');

%%%% Initial analyses - Phase, MAD, BC, Orientation.
% RD-TD correction
rot = rotation.byAxisAngle(zvector,90*degree);

ebsd = rotate(ebsd,rot,'keepXY');
% Index/Phase analysis.
figure;plot(ebsd,'coordinates','on');
ebsdAl=ebsd('Aluminium');
saveas(gcf, 'EBSDmap_ebsd_B_1.png')%%%%%%%%
% BC/MAD analyses
figure; plot(ebsd,ebsd.bc,'coordinates','on');
colormap gray; % make the image grayscale
mtexColorbar;
saveas(gcf, 'EBSDbcmap_ebsd_B_1.png')%%%%%%%%
figure; plot(ebsd,ebsd.mad,'coordinates','on');mtexColorbar;
saveas(gcf, 'EBSDmadmap_ebsd_B_1.png')%%%%%%%%
figure; histogram(ebsd.mad);
saveas(gcf, 'EBSDmaddist_ebsd_B_1.png')%%%%%%%%
figure; histogram(ebsd.bc);
saveas(gcf, 'EBSDbcdist_ebsd_B_1.png')%%%%%%%%

% Plot initial EBSD orientation map.
ipfKey = ipfColorKey(ebsd('Aluminium')); 
ipfKey.inversePoleFigureDirection = vector3d.Z;
figure;plot(ipfKey);% this is the colored fundamental sector
colors = ipfKey.orientation2color(ebsd('Aluminium').orientations);
figure;
plot(ebsd('Aluminium'),colors);
saveas(gcf, 'EBSDOri_ebsd_B_1.png')%%%%%%%%
% Corredction based on MAD distribution 
% Consider only indexed & corrected data.
ebsd_corrected = ebsd(ebsd.mad<1);
ebsdcorrAl=ebsd_corrected('Aluminium');

%% Grain reconstruction.
% Reconstruct the grain structure.
[grains,ebsdcorrAl.grainId,ebsdcorrAl.mis2mean] = calcGrains(ebsdcorrAl,'angle',15*degree);
initialGrainNr=length(grains);
% Delete the very small grains which might be caused by the measurement error.
ebsdcorrAl(grains(grains.grainSize<2)) = []; 
% Redo grain segmentation.
[grains,ebsdcorrAl.grainId] = calcGrains(ebsdcorrAl,'angle',15*degree); 
% Pick up the focused phase.
grainsAl=grains('Aluminium');
totalGrainNr=length(grainsAl);
% Plotting grain mean orientation maps.
figure;
plot(grainsAl,grainsAl.meanOrientation);
saveas(gcf, 'GrianMeanOri_ebsd_B_1.png')%%%%%%%%
% Plotting ebsd orientation map+grain boundaries.
ipfKey = ipfColorKey(ebsdcorrAl); 
ipfKey.inversePoleFigureDirection = vector3d.Z;
colors = ipfKey.orientation2color(ebsdcorrAl.orientations);
figure;plot(ebsdcorrAl,colors);
hold on
plot(grainsAl.boundary,'linewidth',1);
hold off
saveas(gcf, 'EBSDOri+GBs_ebsd_B_1.png')%%%%%%%%

% EBSD KAM map 
% Order 2, distance = 2*step size.
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',2)./degree,'micronbar','off');
mtexColorbar
caxis([0,5])
hold on
plot(grainsAl.boundary,'lineWidth',1.5)
hold off
saveas(gcf, 'GrianKAMmapO2+GBs_ebsd_B_1.png')%%%%%%%%
% Order 3, distance = 3*step size.
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',3)./degree,'micronbar','off');
mtexColorbar
%mtexColorMap LaboTeX
caxis([0,5])
hold on
plot(grainsAl.boundary,'lineWidth',1.5)
hold off
saveas(gcf, 'GrianKAMmapO3+GBs_ebsd_B_1.png')%%%%%%%%

%% Grain reconstruction.
% Reconstruct the grain structure.
[grains,ebsdcorrAl.grainId,ebsdcorrAl.mis2mean] = calcGrains(ebsdcorrAl,'angle',15*degree);
initialGrainNr=length(grains);
% Delete the very small grains which might be caused by the measurement error.
ebsdcorrAl(grains(grains.grainSize<2)) = []; 
% Redo grain segmentation.
[grains,ebsdcorrAl.grainId] = calcGrains(ebsdcorrAl,'angle',15*degree); 
% Pick up the focused phase.
grainsAl=grains('Aluminium');
totalGrainNr=length(grainsAl);
% Plotting grain mean orientation maps.
figure;
plot(grainsAl,grainsAl.meanOrientation);
saveas(gcf, 'GrianMeanOri_ebsd_B_1.png')%%%%%%%%
% Plotting ebsd orientation map+grain boundaries.
ipfKey = ipfColorKey(ebsdcorrAl); 
ipfKey.inversePoleFigureDirection = vector3d.Z;
colors = ipfKey.orientation2color(ebsdcorrAl.orientations);
figure;plot(ebsdcorrAl,colors);
hold on
plot(grainsAl.boundary,'linewidth',1);
hold off
saveas(gcf, 'EBSDOri+GBs_ebsd_B_1.png')%%%%%%%%

% EBSD KAM map 
% Order 2, distance = 2*step size.
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',2)./degree,'micronbar','off');
mtexColorbar
caxis([0,5])
hold on
plot(grainsAl.boundary,'lineWidth',1.5)
hold off
saveas(gcf, 'GrianKAMmapO2+GBs_ebsd_B_1.png')%%%%%%%%
% Order 3, distance = 3*step size.
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',3)./degree,'micronbar','off');
mtexColorbar
%mtexColorMap LaboTeX
caxis([0,5])
hold on
plot(grainsAl.boundary,'lineWidth',1.5)
hold off
saveas(gcf, 'GrianKAMmapO3+GBs_ebsd_B_1.png')%%%%%%%%

%% Grain Size & Shape Data Analysis.
% Find the boundary grains.
outerBoundary_id = any(grainsAl.boundary.grainId==0,2);
grain_id = grainsAl.boundary(outerBoundary_id).grainId;
grain_id(grain_id==0) = [];
% Plot the boundary grains with their mean orientations.
figure;
plot(grainsAl(grain_id),grainsAl(grain_id).meanOrientation);
saveas(gcf, 'BoundaryGrainMap_ebsd_B_1.png')%%%%%%%%
% Remove the boundary grains.
grainsAl(grain_id) = []; 
innerGrainNr=length(grainsAl);
% Plot the inner grains with their mean orientations.
figure;plot(grainsAl,grainsAl.meanOrientation);
saveas(gcf, 'InnerGrainMap_ebsd_B_1.png')%%%%%%%%
% Extract Grain Data.
Grainsize=grainsAl.grainSize;
Grainarea=grainsAl.area;
GraineqR=grainsAl.equivalentRadius;
GraineqD=GraineqR*2;
Grainsf=grainsAl.shapeFactor;
Grainasp=1./grainsAl.aspectRatio;
% Fit the equivalent ellipses of grains.
[GrainfitEangle,GrainfitElongA,GrainfitEshortb] = fitEllipse(grainsAl);
GrainfitEangle=GrainfitEangle./degree;
figure; 
plot(grainsAl,grainsAl.meanOrientation,'linewidth',1);
hold on; 
plotEllipse(grainsAl.centroid,GrainfitElongA,GrainfitEshortb,GrainfitEangle,'lineColor','r');
hold off;
saveas(gcf, 'InnerGrainMap+fitE_ebsd_B_1.png')%%%%%%%%

GrainData=[Grainsize Grainarea GraineqR GraineqD Grainsf GrainfitEangle GrainfitElongA GrainfitEshortb Grainasp];

% Write grain data file
fidw_txt = fopen('GrainData_ebsd_B_1.txt','w+');%%%%%%%%
fprintf(fidw_txt,'There are %6i grains.\r\n',innerGrainNr);
fprintf(fidw_txt,'Grainsize GrainFarea  GrainFeqR GrainFeqD GrainFsf GrainFfitEangle GrainFfitElongA GrainFfitEshortb GrainFasp\r\n');
fprintf(fidw_txt,'%10.2f %12.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\r\n', GrainData');
disp('Grain data writting completed!');
fclose(fidw_txt);

%% Initial grain size/shape relation aanlysis
sz=1;
figure;
scatter(GraineqD,Grainasp,sz);
xlabel('Grain Diameter, {\mu}m','fontsize',15); 
ylabel('Grain Shape Aspect Ratio, - ','fontsize',15);
saveas(gcf, 'GrainSize+ShapeAsp_ebsd_B_1.png')%%%%%%%%
figure;
scatter(GraineqD,GrainfitEangle,sz);
xlabel('Grain Diameter, {\mu}m','fontsize',15); 
ylabel('Grain Shape Theta angle, degree','fontsize',15);
ylim([0 180]);
saveas(gcf, 'GrainSize+ShapeAngle_ebsd_B_1.png')%%%%%%%%
figure;
scatter(Grainasp,GrainfitEangle,sz);
xlabel('Grain Shape Aspect Ratio, - ','fontsize',15); 
ylabel('Grain Shape Theta angle, degree','fontsize',15);
ylim([0 180]);
saveas(gcf, 'GrainShapeAsp+Angle_ebsd_B_1.png')%%%%%%%%

%% Grain size - number fraction distribution
% Plot data originally in dataset "diameter data".
[CdfF,CdfX] = ecdf(GraineqD,'Function','cdf');% compute empirical cdf
BinInfo.rule = 5;
BinInfo.width = 1;
BinInfo.placementRule = 1;
[~,BinEdge] = internal.stats.histbins(GraineqD,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
figure;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','Grain Size Data','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1);
xlabel('Grain Diameter, {\mu}m','fontsize',15); 
ylabel('Number Fraction,-','fontsize',15);
% Create grid where function will be computed.
hold on;
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
% Fit this distribution to get parameter values.
pdnum = fitdist(GraineqD,'lognormal');
mu_n=pdnum.mu;
sigma_n=pdnum.sigma;
mean_n = mean(pdnum);
median_n = median(pdnum);
v_n = std(pdnum);
mode_n=exp(mu_n-sigma_n.^2);
YPlot = pdf(pdnum,XGrid);
hLine = plot(XGrid,YPlot,'Color','b','LineStyle','-', 'LineWidth',1.5,'DisplayName','Log-normal distribution');
legend('show','Location','northeast')
hold off;
saveas(gcf, 'GrainSizeNumberLog_ebsd_B_1.png')%%%%%%%%


%% Grain size - area fraction distribution
% Plot data originally in dataset "diameter data" + grianSize as weights
weight=grainsAl.grainSize;
[CdfF,CdfX] = ecdf(GraineqD,'Function','cdf','freq',weight);  % compute empirical cdf
BinInfo.width = 1;
[~,BinEdge] = internal.stats.histbins(GraineqD,[],weight,BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
figure;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','Grain Size Data','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1);
xlabel('Grain Diameter, {\mu}','fontsize',15); 
ylabel('Area Fraction,-','fontsize',15);
hold on;
% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);

% Fit this distribution to get parameter values
pdarea = fitdist(GraineqD,'lognormal','freq',weight);
mu_a=pdarea.mu;
sigma_a=pdarea.sigma;
mean_a = mean(pdarea);
median_a = median(pdarea);
v_a = std(pdarea);
mode_a=exp(mu_a-sigma_a.^2);
YPlot = pdf(pdarea,XGrid);
hLine = plot(XGrid,YPlot,'Color','b','LineStyle','-', 'LineWidth',1.5,'DisplayName','Log-normal distribution');
legend('show','Location','northeast')
hold off;
saveas(gcf, 'GrainSizeAreaLog_ebsd_B_1.png')%%%%%%%%

%% Grain shape aspect ratio - number fraction distribution
% Plot data originally in dataset "shape data".
[CdfF,CdfX] = ecdf(Grainasp,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 5;
BinInfo.width = 0.02;
BinInfo.placementRule = 1;
[~,BinEdge] = internal.stats.histbins(Grainasp,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
BinHeight=BinHeight*0.05;
figure;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','Grain Shape Data','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1);
xlabel('Grain Shape Aspect Ratio, - ','fontsize',15);
ylabel('Number Fraction, -','fontsize',15)
hold on
% Create grid where function will be computed.
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
% Fit this distribution to get parameter values - Beta.
pdaspBeta = fitdist(Grainasp,'Beta');
a_aspBeta = pdaspBeta.a;
b_aspBeta = pdaspBeta.b;
mean_aspBeta = mean(pdaspBeta);
v_aspBeta = std(pdaspBeta);
median_aspBeta = median(pdaspBeta);
YPlot = pdf(pdaspBeta,XGrid)*0.05;
hLine = plot(XGrid,YPlot,'Color','r','LineStyle','-', 'LineWidth',1.5,'DisplayName','Beta distribution');
legend('show','Location','northwest')
xlim([0 1]);
hold off;
saveas(gcf, 'GrainAspNumberBeta_ebsd_B_1.png')%%%%%%%%

%% Grain shape angle - number fraction distribution
[CdfF,CdfX] = ecdf(GrainfitEangle,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 5;
BinInfo.width = 5;
BinInfo.placementRule = 1;
[~,BinEdge] = internal.stats.histbins(GrainfitEangle,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
BinHeight=BinHeight*0.05;
figure;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','Grain Shape Data','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1);
xlabel('Grain Shape Theta angle, degree ','fontsize',15);
xlim([0 180]);
xticks(0:30:180);
ylabel('Number Fraction, -','fontsize',15)
legend('show','Location','northwest')
saveas(gcf, 'GrainShapeAngleNumber_ebsd_B_1.png')%%%%%%%%


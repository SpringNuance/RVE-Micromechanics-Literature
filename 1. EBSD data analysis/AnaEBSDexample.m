%% Script for EBSD Data analysis
% Purpose: EBSD analysis with Matlab-MTEX toolbox
% Functions: Import EBSD data
%            EBSD plotting & initial analyses
%            Grain reconstruction
%            Grain Size & Shape analysis
%            Grain Orientation analysis
%            Grain Misorientation analysis
%
% Record of revision:
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  30-10-2020        Wenqi Liu       Original code for MTEX 5.4.0
%
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  26-11-2020        Wenqi Liu        Modify for Al measurement
%
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  26-11-2020        Wenqi Liu        Modify MisOri codes, 
%                                   Check: problem with MisOri2!
%
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  27-11-2020        Wenqi Liu        Change X/Y for RD/TD!!! 
%
%----------------------- General Comments END---------------------------
%% Import the Data.	
% Specify crystal symmetry.
CS = {... 
  'notIndexed',...
  crystalSymmetry('6/mmm', [5.8 5.8 4.6], 'X||a*', 'Y||b', 'Z||c', 'mineral', 'Ti3Al - alpha2', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('4/mmm', [4 4 4.1], 'mineral', 'TiAl - gamma', 'color', [0.56 0.74 0.56]),...
  crystalSymmetry('m-3m', [3.3 3.3 3.3], 'mineral', 'Titanium-Cubic', 'color', [0.85 0.65 0.13]),...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.94 0.5 0.5]),...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0 0 0.55])};
% Plotting convention.
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');
% Path to files.
pname = 'C:\Users\liuw7\OneDrive\15 Sheila\100 reference\EBSDDataAna';
% Which files to be imported.
fname = [pname '\RD_TD1.ctf'];
% Create an EBSD variable containing the data.
ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
 
%% Initial analyses - Phase, MAD, BC, Orientation.
% RD-TD correction
rot = rotation('Euler',90*degree,0*degree,0*degree);
ebsd = rotate(ebsd,rot,'keepXY');
% Index/Phase analysis.
figure;plot(ebsd,'coordinates','on');
ebsdAl=ebsd('Aluminium');
saveas(gcf, 'EBSDmap_ebsdRDTD_1.png')
% BC/MAD analyses
figure; plot(ebsd,ebsd.bc,'coordinates','on');
colormap gray; % make the image grayscale
mtexColorbar;
saveas(gcf, 'EBSDbcmap_ebsdRDTD_1.png')
figure; plot(ebsd,ebsd.mad,'coordinates','on');mtexColorbar;
saveas(gcf, 'EBSDmadmap_ebsdRDTD_1.png')
figure; histogram(ebsd.mad);
saveas(gcf, 'EBSDmaddist_ebsdRDTD_1.png')
figure; histogram(ebsd.bc);
saveas(gcf, 'EBSDbcdist_ebsdRDTD_1.png')

% Plot initial EBSD orientation map.
ipfKey = ipfColorKey(ebsd('Aluminium')); 
ipfKey.inversePoleFigureDirection = vector3d.Z;
figure;plot(ipfKey);% this is the colored fundamental sector
colors = ipfKey.orientation2color(ebsd('Aluminium').orientations);
figure;
plot(ebsd('Aluminium'),colors);
saveas(gcf, 'EBSDOri_ebsdRDTD_1.png')
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
saveas(gcf, 'GrianMeanOri_ebsdRDTD_1.png')
% Plotting ebsd orientation map+grain boundaries.
ipfKey = ipfColorKey(ebsdcorrAl); 
ipfKey.inversePoleFigureDirection = vector3d.Z;
colors = ipfKey.orientation2color(ebsdcorrAl.orientations);
figure;plot(ebsdcorrAl,colors);
hold on
plot(grainsAl.boundary,'linewidth',1);
hold off
saveas(gcf, 'EBSDOri+GBs_ebsdRDTD_1.png')

% EBSD KAM map
% Order 2, distance = 2*step size.
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',2)./degree,'micronbar','off');
mtexColorbar
caxis([0,5])
hold on
plot(grainsAl.boundary,'lineWidth',1.5)
hold off
saveas(gcf, 'GrianKAMmapO2+GBs_ebsdRDTD_1.png')
% Order 3, distance = 3*step size.
figure;plot(ebsd,ebsd.KAM('threshold',5*degree,'order',3)./degree,'micronbar','off');
mtexColorbar
%mtexColorMap LaboTeX
caxis([0,5])
hold on
plot(grainsAl.boundary,'lineWidth',1.5)
hold off
saveas(gcf, 'GrianKAMmapO3+GBs_ebsdRDTD_1.png')

%% Grain Size & Shape Data Analysis.
% Find the boundary grains.
outerBoundary_id = any(grainsAl.boundary.grainId==0,2);
grain_id = grainsAl.boundary(outerBoundary_id).grainId;
grain_id(grain_id==0) = [];
% Plot the boundary grains with their mean orientations.
figure;
plot(grainsAl(grain_id),grainsAl(grain_id).meanOrientation);
saveas(gcf, 'BoundaryGrainMap_ebsdRDTD_1.png')
% Remove the boundary grains.
grainsAl(grain_id) = []; 
innerGrainNr=length(grainsAl);
% Plot the inner grains with their mean orientations.
figure;plot(grainsAl,grainsAl.meanOrientation);
saveas(gcf, 'InnerGrainMap_ebsdRDTD_1.png')
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
saveas(gcf, 'InnerGrainMap+fitE_ebsdRDTD_1.png') 

GrainData=[Grainsize Grainarea GraineqR GraineqD Grainsf GrainfitEangle GrainfitElongA GrainfitEshortb Grainasp];

% Write grain data file
fidw_txt = fopen('GrainData_ebsdRDTD_1.txt','w+');
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
saveas(gcf, 'GrainSize+ShapeAsp_ebsdRDTD_1.png')
figure;
scatter(GraineqD,GrainfitEangle,sz);
xlabel('Grain Diameter, {\mu}m','fontsize',15); 
ylabel('Grain Shape Theta angle, degree','fontsize',15);
ylim([0 180]);
saveas(gcf, 'GrainSize+ShapeAngle_ebsdRDTD_1.png')
figure;
scatter(Grainasp,GrainfitEangle,sz);
xlabel('Grain Shape Aspect Ratio, - ','fontsize',15); 
ylabel('Grain Shape Theta angle, degree','fontsize',15);
ylim([0 180]);
saveas(gcf, 'GrainShapeAsp+Angle_ebsdRDTD_1.png')

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
saveas(gcf, 'GrainSizeNumberLog_ebsdRDTD_1.png')


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
saveas(gcf, 'GrainSizeAreaLog_ebsdRDTD_1.png')

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
saveas(gcf, 'GrainAspNumberBeta_ebsdRDTD_1.png')

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
saveas(gcf, 'GrainShapeAngleNumber_ebsdRDTD_1.png')

%% Write parameters
fidw_txt = fopen('GrainDataParameters_ebsdRDTD_1.txt','w+');
fprintf(fidw_txt,'AL RD-TD-1 plane\r\n');
fprintf(fidw_txt,'Lognormal distribution:\r\n');
fprintf(fidw_txt,'       mu             sigma          mean           median         mode           sqrt\r\n');
fprintf(fidw_txt,'d-num  %12.8f%15.8f%15.8f%15.8f%15.8f%15.8f\r\n',[mu_n sigma_n mean_n median_n mode_n v_n]);
fprintf(fidw_txt,'d-area %12.8f%15.8f%15.8f%15.8f%15.8f%15.8f\r\n',[mu_a sigma_a mean_a median_a mode_a v_a]);

fprintf(fidw_txt,'\r\n\r\nBeta distribution:\r\n');
fprintf(fidw_txt,'       a              b              mean           median         sqrt\r\n');
fprintf(fidw_txt,'asp-num   %12.8f%15.8f%15.8f%15.8f%15.8f\r\n',[a_aspBeta b_aspBeta mean_aspBeta median_aspBeta v_aspBeta]);
disp('Grain data parameter writting completed!');
fclose(fidw_txt);

%% %% Grain Orientation & misorientation analyses
% Based on graindata-del2
psidelS=calcKernel(grainsAl.meanOrientation);
odfdelS=calcDensity(ebsdcorrAl.orientations,'kernel',psidelS,'Bunge','ZXZ');
odfdelS.SS=specimenSymmetry('-1');
h=[Miller(1,0,0,odfdelS.CS),Miller(1,1,0,odfdelS.CS),Miller(1,1,1,odfdelS.CS)];
tdelS=textureindex(odfdelS,'resolution',5*degree); 
% Pole figures plotting
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'RD','TD'},...
 'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','intoPlane');
figure; 
plotPDF(odfdelS,h,'contourf','resolution',8*degree,'MinMax'); 
CLim(gcm,'equal');setColorRange([0 3],'current');mtexColorbar;
saveas(gcf, 'ebsdRDTD_1_m-3m_PF_GB15del2.png')
% Inverse pole figures plotting
figure;plotIPDF(odfdelS,[xvector,yvector,zvector],'antipodal','MinMax');
CLim(gcm,'equal');setColorRange([0 2],'current');mtexColorbar; 
saveas(gcf, 'ebsdRDTD_1_m-3m_IPF_GB15del2.png')
% ODF figures plotting
figure; 
odfdelS.SS=specimenSymmetry('orthorhombic');
plot3d(odfdelS,'MinMax');mtexColorbar;
saveas(gcf, 'ebsdRDTD_1_m-3m_3DODF_GB15del2.png')
figure;plot(odfdelS,'resolution',8*degree,'sections',18,'projection','plain','contourf','MinMax');
CLim(gcm,'equal');setColorRange([0 4],'current');mtexColorbar;
saveas(gcf, 'ebsdRDTD_1_m-3m_ODF_GB15del2.png')
% Grain misorientation (uncorrelated)
mori1= calcMisorientation(ebsdcorrAl,ebsdcorrAl);
misorientation1=angle(mori1)./degree;
[n1,central1] = hist(misorientation1,15:2.39:62.8);
morifraction1 = n1./sum(n1);
figure;
plot(central1,morifraction1,'b','LineWidth',1);
xlabel('Misorientation,- ','fontsize',15);
ylabel('Number Fraction,-','fontsize',15);
legend('uncorrelatedMisO','Location','northwest','fontsize',15);
saveas(gcf, 'ebsdRDTD_1_m-3m_misOriDis1_GB15del2.png')

% % Grain misorientation (correlated)
% GBsgrainscorri=grains.boundary;
% % GBsAL=GBsgrainscorri('Aluminium','Aluminium');
% % misorientation2 = GBsAL.misorientation.angle./degree;
% rots=rotation.byEuler(GBsgrainscorri.misrotation.phi1,GBsgrainscorri.misrotation.Phi,GBsgrainscorri.misrotation.phi2);
% misorientation2 = angle(rots)./degree;
% [n2,central2] = hist(misorientation2,15:2.39:62.8);
% morifraction2 = n2./sum(n2);
% figure;
% plot(central1,morifraction1,'b','LineWidth',1);
% hold on; 
% plot(central2,morifraction2,'r--','LineWidth',1);
% hold off;
% xlim([0 62.8]);
% xlabel('Misorientation,- ','fontsize',15);
% ylabel('Number Fraction,-','fontsize',15);
% legend({'uncorrelatedMisO','correlatedMisO'},'Location','northwest','fontsize',15);
% saveas(gcf, 'ebsdRDTD_1_m-3m_misOriDis_GB15del2.png')
% figure;
% plot(GBsgrainscorri,misorientation2,'linewidth',3)
% mtexColorbar('title','misorientation angle (°)')

%% Wirte EBSD orientation data file for Dream3D input
ebsdOris=ebsdcorrAl.orientations;
Nori=length(ebsdOris);
EulerAng=ones(Nori,5);
EulerAng(:,1:3)=[ebsdOris.phi1 ebsdOris.Phi ebsdOris.phi2];
fidw_txt = fopen('ebsdRDTD_1_corriAldel2_Oris_RVE.txt','w+');
fprintf(fidw_txt,'# Euler0 Euler1 Euler2 Weight Sigma\r\n');
fprintf(fidw_txt,'Angle Count:%i\r\n',Nori);
fprintf(fidw_txt,'%7.5f%9.5f%9.5f%3.i%3.i\r\n',EulerAng');
disp('Orientation file for RVE writting completed!');
fclose(fidw_txt);
function [] = DREAME3DRVEoptGS(ElementNumber,MeshSize,RVEID)
% Script file: DREAME3DRVEoptGS.m
%
% Purpose: 1. Analyze RVE Grain size distribution.
%          2. Analyze RVE Grain Shape distribution.
%
% Note: 1. Import txt file name shall be: RVEFFT_a_b_c.txt 
%			a:ElementNumber b:MeshSize c:RVEID, e.g. RVEFFT_40_2_1.txt
%		2. Import csv file name shall be: RVE_a_b_c.csv 
%			a:ElementNumber b:MeshSize c:RVEID, e.g. RVE_40_2_1.csv
%		3. The Exp./Input grain SIZE/SHAPE parameters shall be given in lines 45-50.
%
% Record of rivisions:
%     Date           Programmer        Description of change 
%   ========         ==========     ===========================
%  09-07-2017        Wenqi Liu           Original code
%
%     Date           Programmer        Description of change 
%   ========         ==========     ===========================
%  11-03-2020        Wenqi Liu       1.Delete other CP inputs 
%                                    2.Indicate the length unit
%
%     Date           Programmer        Description of change 
%   ========         ==========     ===========================
%  02-04-2020        Wenqi Liu      1.Add the grain size analysis
%                                   2.Correct the mesh size to box size
%                                   3.Add grain size comparison with
%                                   Exp.data
%
%     Date           Programmer        Description of change 
%   ========         ==========     ===========================
%  03-07-2020        Wenqi Liu           Change for nm
%
%     Date           Programmer        Description of change 
%   ========         ==========     ===========================
%  11-12-2020        Wenqi Liu      Change input as um & output to mm
%
%     Date           Programmer        Description of change 
%   ========         ==========     ===========================
%  21-12-2020        Wenqi Liu      1. Delete input writting
%									2. Add grian shape mean value comparison
%----------------------- General Comments END---------------------------
%
%% Input data from Exp. 
expMuN = 6.2409727;    %%%%% Exp. NF lognormal distribution mu 
expSigmaN = 0.69440648;%%%%% Exp. NF lognormal distribution sigma 
expMuA = 7.60216176;   %%%%% Exp. AF lognormal distribution mu 
expSigmaA = 0.78488038;%%%%% Exp. NF lognormal distribution sigma 
AspExpRDTD=0.59728906;      %%%%%%%%%%%%%Exp. Shape Asp
AspExpRDND=0.59728906;      %%%%%%%%%%%%%Exp. Shape Asp

%% Initialize variables.
nElement_x = ElementNumber;   %%%% element number_x
nElement_y = ElementNumber;   %%%% element number_y
nElement_z = ElementNumber;   %%%% element number_z
% MESH SIZE UNITE: nm!
mesh_x = MeshSize;   %%%% gridSize_x in um!
mesh_y = MeshSize;   %%%% gridSize_y in um!
mesh_z = MeshSize;   %%%% gridSize_z in um!
mesh_x_mm = mesh_x/1000; %%%% change RVE unit: um2mm 
mesh_y_mm = mesh_y/1000; %%%% change RVE unit: um2mm 
mesh_z_mm = mesh_z/1000; %%%% change RVE unit: um2mm 
% BOX/RVE SIZE UNITE: mm!
Size_x = nElement_x*mesh_x_mm;   %%%% RVESize_x in mm
Size_y = nElement_y*mesh_y_mm;   %%%% RVESize_y in mm
Size_z = nElement_z*mesh_z_mm;   %%%% RVESize_z in mm

%% Import data from FFT.txt
filename = sprintf('RVEFFT_%d_%d_%d.txt',[ElementNumber MeshSize RVEID]);%%%%% File name
delimiter = ' ';
% Format string for each line of text:
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of FFT_Dream3d according to format string.
FFT_Dream3dArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Create output variable
FFT_Dream3d = [FFT_Dream3dArray{1:end-1}];
%%%%% [Phi1] [PHI] [Phi2] [X] [Y] [Z] [Grain Id] [Phase Id]
% Clear temporary variables
clearvars filename delimiter formatSpec fileID FFT_Dream3dArray ans;

%% Create the grain based FFT_Dream3d
nGrain =  max(FFT_Dream3d(:,7));
F = cell(nGrain,1);
grains = zeros(nGrain,8);
grains(:,1) = 1:1:nGrain;          %%%%Grain ID
for i = 1:1:nGrain
    grains(i,6) = length(find(FFT_Dream3d(:,7)==i));
    F{i,1} = find(FFT_Dream3d(:,7)==i);
    grains(i,2) = FFT_Dream3d((F{i,1}(1,1)),8);     %%%%Grain: Phase ID
    grains(i,3) = FFT_Dream3d((F{i,1}(1,1)),1);     %%%%Grain: Euler phi1 in degree
    grains(i,4) = FFT_Dream3d((F{i,1}(1,1)),2);     %%%%Grain: Euler Phi in degree
    grains(i,5) = FFT_Dream3d((F{i,1}(1,1)),3);     %%%%Grain: Euler phi2 in degree
end
grains(:,7) = grains(:,6)/sum(grains(:,6));  % Grain volumn fraction

% % If RVE has two phases.
% Phase1=grains(find(grains(:,2)==1),8);
% Phase2=grains(find(grains(:,2)==2),8);

vElement=mesh_x*mesh_y*mesh_z;
Pgrain=grains(:,6); % Grain points number
Vgrain=grains(:,6)*vElement; % Grain volumn
Rgrain=(3*Vgrain/(4*pi)).^(1/3); % Grain radius
Dgrain=2*Rgrain;    % Grain diameter in um
% Dgrain_mm=Dgrain/1000;    % Grain diameter in mm
grains(:,8)=Dgrain; % Grain diameter
Agrain=Rgrain.^2*pi; % Grain area

% Write grain size data in RVE
Gfilename=sprintf('RVE%d_%d_%d_GrainSizeData.txt',[ElementNumber MeshSize RVEID]);
fidw_txt = fopen(Gfilename,'w+');
fprintf(fidw_txt,'RVE%d_%d_%d_GrainSizeData\r\n',[ElementNumber MeshSize RVEID]);
fprintf(fidw_txt,'Diameter    Volumn     PointsN\r\n');
fprintf(fidw_txt,'%15.8f%15.8f%12i\r\n',[Dgrain Vgrain Pgrain]');
disp('Grain size data writting completed!');
fclose(fidw_txt);

%% Grain size distribution analysis
Dmax=ceil(max(Dgrain)/10)*10; % Set maximum diamater

% Number fraction
% Create grid where function will be computed
[CdfF,CdfX] = ecdf(Dgrain,'Function','cdf');% compute empirical cdf
BinInfo.rule = 5;
BinInfo.width = 200;
BinInfo.placementRule = 1;
[~,BinEdge] = internal.stats.histbins(Dgrain,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
figure;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','RVE GD:NF','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1.5);
xlabel('Grain diameter, um','fontsize',15); 
ylabel('Number fraction,-','fontsize',15);
xlim([0 Dmax]);
hold on;
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
% Number fraction fitting
pdnum = fitdist(Dgrain,'lognormal');
mu_n=pdnum.mu;
sigma_n=pdnum.sigma;
mean_n = mean(pdnum);
median_n = median(pdnum);
v_n = std(pdnum);
mode_n=exp(mu_n-sigma_n.^2);
YPlot = pdf(pdnum,XGrid);
hLine = plot(XGrid,YPlot,'Color','b','LineStyle','-', 'LineWidth',1.5,'DisplayName','Log-normal:NF');
legend('show','Location','northeast')
hold off;
figurename = sprintf('RVE%d_%d_%d_GSNF.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
% Volumn fraction
% Create grid where function will be computed
[CdfF,CdfX] = ecdf(Dgrain,'Function','cdf','freq',Pgrain);  % compute empirical cdf
BinInfo.width = 200;
[~,BinEdge] = internal.stats.histbins(Dgrain,[],Pgrain,BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
figure;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','RVE GD:VF','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1.5);
xlabel('Grain diameter, um','fontsize',15); 
ylabel('Volumn fracion,-','fontsize',15);
xlim([0 Dmax]);
hold on;
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
% Volumn fraction fitting
pdvolumn = fitdist(Dgrain,'lognormal','freq',Pgrain);
mu_v = pdvolumn.mu;
sigma_v = pdvolumn.sigma;
mean_v = mean(pdvolumn);
median_v = median(pdvolumn);
v_v = std(pdvolumn);
mode_v=exp(mu_v-sigma_v.^2);
YPlot = pdf(pdvolumn,XGrid);
hLine = plot(XGrid,YPlot,'Color','b','LineStyle','-', 'LineWidth',1.5,'DisplayName','Log-normal:VF');
legend('show','Location','northwest')
hold off;
figurename = sprintf('RVE%d_%d_%d_GSVF.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
GrainSize=[mu_n sigma_n mean_n median_n mode_n v_n;...
          mu_v sigma_v mean_v median_v mode_v v_v];

Gfilename=sprintf('RVE_GrainSizeParameters_um.txt');
fidw_txt = fopen(Gfilename,'a+');
fprintf(fidw_txt,'RVE%d_%d_%d\r\n',[ElementNumber MeshSize RVEID]);
fprintf(fidw_txt,'Lognormal distribution:\r\n');
fprintf(fidw_txt,'       mu             sigma          mean           median         mode           sqrt\r\n');
fprintf(fidw_txt,'d-num  %15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\r\n',[mu_n sigma_n mean_n median_n mode_n v_n]);
fprintf(fidw_txt,'d-vol  %15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\r\n',[mu_v sigma_v mean_v median_v mode_v v_v]);
fprintf(fidw_txt,'\r\n');
disp('Grain size parameter writting completed!');
fclose(fidw_txt);
%% Comparison of grain size distribution 
xdia = 0:(Dmax/100):Dmax;
fitExpN=lognpdf(xdia,expMuN,expSigmaN);
fitExpA=lognpdf(xdia,expMuA,expSigmaA);
fitRVEN=lognpdf(xdia,mu_n,sigma_n);
fitRVEV=lognpdf(xdia,mu_v,sigma_v);
% Fitting curves
figure; hold on;
plot(xdia,fitExpN,'k','DisplayName','Exp.NF','linewidth',1.5);
plot(xdia,fitExpA,'k--','DisplayName','Exp.AF','linewidth',1.5);
plot(xdia,fitRVEN,'b','DisplayName','RVE.NF','linewidth',1.5);
plot(xdia,fitRVEV,'b--','DisplayName','RVE.VF','linewidth',1.5);
xlabel('Grain diameter, um','fontsize',15); ylabel('Fraction,-','fontsize',15);
legend('Location','northeast');box('on');
set(gca,'fontsize',15,'fontname','Arial');
hold off;
figurename = sprintf('RVE%d_%d_%d_GS3.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
% Exp. fitting curves + RVE data
[CdfF,CdfX] = ecdf(Dgrain,'Function','cdf');% compute empirical cdf
BinInfo.rule = 5;
BinInfo.width = 200;
BinInfo.placementRule = 1;
[~,BinEdge] = internal.stats.histbins(Dgrain,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
figure;hold on;
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','RVE GD:NF','FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',1.5);

[CdfF,CdfX] = ecdf(Dgrain,'Function','cdf','freq',Pgrain);  % compute empirical cdf
BinInfo.width = 200;
[~,BinEdge] = internal.stats.histbins(Dgrain,[],Pgrain,BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'DisplayName','RVE GD:VF','FaceColor','none','EdgeColor','b','LineStyle','-', 'LineWidth',1.5);
xlabel('Grain diameter, um','fontsize',15); 
ylabel('Fraction,-','fontsize',15);
xlim([0 Dmax]);

XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);
fitExpN=lognpdf(XGrid,expMuN,expSigmaN);
fitExpA=lognpdf(XGrid,expMuA,expSigmaA);
plot(XGrid,fitExpN,'k','DisplayName','Exp.NF','linewidth',1.5);
plot(XGrid,fitExpA,'b','DisplayName','Exp.AF','linewidth',1.5);
xlabel('Grain diameter, um','fontsize',15); ylabel('Fraction,-','fontsize',15);
legend('Location','northeast');box('on');
set(gca,'fontsize',15,'fontname','Arial');
hold off;
figurename = sprintf('RVE%d_%d_%d_GS4.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
disp('Grain size distribution plotting completed!');

%% Grain shape mean value comparison
filename = sprintf('RVE_%d_%d_%d.csv',[ElementNumber MeshSize RVEID]);
RVEAsp=readmatrix(filename);
AvgAsp12=[mean(RVEAsp(:,6)) mean(RVEAsp(:,7))];
AvgAspall=mean(AvgAsp12)
AspP=[AvgAsp12 AvgAspall];
Afilename=sprintf('RVE_GrainShapeAspData.txt');
fidw_txt = fopen(Afilename,'a+');
fprintf(fidw_txt,'For material RD-TD plane\r\n');
fprintf(fidw_txt,'Initial Grain Asp = %12.8f\r\n',AspExpRDTD);
fprintf(fidw_txt,'\r\n');
fprintf(fidw_txt,'For material RD-ND plane\r\n');
fprintf(fidw_txt,'Initial Grain Asp = %12.8f\r\n',AspExpRDND);
fprintf(fidw_txt,'\r\n');
fprintf(fidw_txt,'RVE_%d_%d_%d_GrainAsp\r\n',[ElementNumber MeshSize RVEID]);
fprintf(fidw_txt,' RVEID AvgAspAxis1 AvgAspAxis2 AvgAspAll\r\n');
fprintf(fidw_txt,'%5i%12.8f%12.8f%12.8f\r\n',[RVEID AspP]);
fclose(fidw_txt);
txt=sprintf('RVE_GrainShapeAspData writting completed!');
disp(txt);
end


% Script file: RVEFFToriPlotting.m
%
% Purpose: Analyze scatted orientations FOR SINGLE PHASE!
% Functions: Import orientation data from Dream3D FFT.txt
%            Orientation analyses
%
% Note: based on the inital ebsd analyses as inputs.
%
% Record of rivisions:
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  01-04-2020        Wenqi Liu       Original code for MTEX 5.2.8
%
%     Date           Programmer          Description of change 
%   ========         ==========     =================================
%  03-07-2020        Wenqi Liu          Change import formate
%                                       Add saving figures
%----------------------- General Comments END---------------------------

%% Import data from FFT.txt
ElementNumber=40;
MeshSize=200;
RVEID=1;
filename = sprintf('RVEFFT_%d_%d_%d.txt',[ElementNumber MeshSize RVEID]);%%%%% File name:RVEFFT_40_05_1
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

a=FFT_Dream3d(:,1:3);

%% Construction of ODF
% Defination of crystal and sample symmetry
cs = CS(2);  %% Use the same CS as initial ebsd setting
ss = specimenSymmetry('-1');
% Converting data to orientation
Ori=orientation('Euler',a*degree,odfini.CS); %% Use the same odfini as initial ebsd setting
%rot=rotation('euler',0*degree,90*degree,0*degree);%%-------if necessary
%ori2=rot*ori1;
% Construction of ODF
odf=calcODF(Ori,'kernel',psiini,'Bunge','ZXZ','silent'); %% Use the same psiini as initial ebsd setting
t=textureindex(odf,'resolution',5*degree);
e=calcError(odf,odfini,'L2'); %% calculate error between initial ODF and RVE ODF
ans=[t e]

%% Figures plotting
% Plot Pole Figures
h=[Miller(1,0,0,odfini.CS),Miller(1,1,0,odfini.CS),Miller(1,1,1,odfini.CS)];
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','intoPlane');
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'RD','TD'},...
  'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
odf.SS=specimenSymmetry('-1');
figure;plotPDF(odf,h,'contourf','resolution',8*degree,'MinMax'); 
CLim(gcm,'equal');setColorRange([0 3],'current');mtexColorbar;
figurename = sprintf('PF_GrainDelS5GB15_RVE%d_%d_%d.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
% Plot Inverse Pole Figures
figure;plotIPDF(odf,[xvector,yvector,zvector],'antipodal','MinMax');
CLim(gcm,'equal');setColorRange([0 2],'current');mtexColorbar; 
figurename = sprintf('IPF_GrainDelS5GB15_RVE%d_%d_%d.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
% Plot ODF Figures
odf.SS=specimenSymmetry('orthorhombic');
figure;plot(odf,'resolution',8*degree,'sections',18,'projection','plain',...
    'contourf','MinMax','silent');
CLim(gcm,'equal');setColorRange([0 7],'current');mtexColorbar;
figurename = sprintf('ODF_GrainDelS5GB15_RVE%d_%d_%d.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)
figure;plot3d(odf);mtexColorbar;
figurename = sprintf('3DODF_GrainDelS5GB15_RVE%d_%d_%d.png',[ElementNumber MeshSize RVEID]);
saveas(gcf, figurename)


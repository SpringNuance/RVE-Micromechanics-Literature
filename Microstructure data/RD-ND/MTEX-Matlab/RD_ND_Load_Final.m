%% Import Script for EBSD Data

clear,clc

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('6/mmm', [5.8 5.8 4.6], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti3Al - alpha2', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('4/mmm', [4 4 4.1], 'mineral', 'TiAl - gamma', 'color', [0.56 0.74 0.56]),...
  crystalSymmetry('m-3m', [3.3 3.3 3.3], 'mineral', 'Titanium-Cubic', 'color', [0.85 0.65 0.13]),...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.94 0.5 0.5]),...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0 0 0.55])};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('yAxisDirection','north');


%% IMPORT THE DATA

%% SAMPLE 1

% create an EBSD variable containing the data
ebsd1 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND1.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
P1=sum(ebsd1.isIndexed == 1);
odf1 = calcODF(ebsd1('aluminium').orientations);

grains1 = calcGrains(ebsd1('indexed'),'threshold',15*degree)
grains1 = smooth(grains1,5);
% Create border grain array
i=1;
Border1=[];
while i<=numel(grains1.boundary.grainId)/2
    if grains1.boundary.grainId(i,1)==0
        Border1(end+1)=grains1.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains1 = grains1;
selected_grains1(Border1)=[]
[omega1,a1,b1] = selected_grains1.fitEllipse;

%% SAMPLE 2

% create an EBSD variable containing the data
ebsd2 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND2.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd2.prop.y=ebsd2.prop.y+402;
P2=sum(ebsd2.isIndexed == 1);
odf2 = calcODF(ebsd2('aluminium').orientations);

grains2 = calcGrains(ebsd2('indexed'),'threshold',15*degree)
grains2 = smooth(grains2,5);
% Create border grain array
i=1;
Border2=[];
while i<=numel(grains2.boundary.grainId)/2
    if grains2.boundary.grainId(i,1)==0
        Border2(end+1)=grains2.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains2 = grains2;
selected_grains2(Border2)=[]
[omega2,a2,b2] = selected_grains2.fitEllipse;


%% SAMPLE 3

% create an EBSD variable containing the data
ebsd3 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND3.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd3.prop.y=ebsd3.prop.y+402+279;
P3=sum(ebsd3.isIndexed == 1);
odf3 = calcODF(ebsd3('aluminium').orientations);

grains3 = calcGrains(ebsd3('indexed'),'threshold',15*degree)
grains3 = smooth(grains3,5);
% Create border grain array
i=1;
Border3=[];
while i<=numel(grains3.boundary.grainId)/2
    if grains3.boundary.grainId(i,1)==0
        Border3(end+1)=grains3.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains3 = grains3;
selected_grains3(Border3)=[]
[omega3,a3,b3] = selected_grains3.fitEllipse;

%% SAMPLE 4

% create an EBSD variable containing the data
ebsd4 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND4.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd4.prop.y=ebsd4.prop.y+402+279*2;
P4=sum(ebsd4.isIndexed == 1);
odf4 = calcODF(ebsd4('aluminium').orientations);

grains4 = calcGrains(ebsd4('indexed'),'threshold',15*degree)
grains4 = smooth(grains4,5);
% Create border grain array
i=1;
Border4=[];
while i<=numel(grains4.boundary.grainId)/2
    if grains4.boundary.grainId(i,1)==0
        Border4(end+1)=grains4.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains4 = grains4;
selected_grains4(Border4)=[]
[omega4,a4,b4] = selected_grains4.fitEllipse;


%% SAMPLE 5

% create an EBSD variable containing the data
ebsd5 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND5.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd5.prop.y=ebsd5.prop.y+402+279*3;
P5=sum(ebsd5.isIndexed == 1);
odf5 = calcODF(ebsd5('aluminium').orientations);

grains5 = calcGrains(ebsd5('indexed'),'threshold',15*degree)
grains5 = smooth(grains5,5);
% Create border grain array
i=1;
Border5=[];
while i<=numel(grains5.boundary.grainId)/2
    if grains5.boundary.grainId(i,1)==0
        Border5(end+1)=grains5.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains5 = grains5;
selected_grains5(Border5)=[]
[omega5,a5,b5] = selected_grains5.fitEllipse;

%% SAMPLE 6

% create an EBSD variable containing the data
ebsd6 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND6.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd6.prop.y=ebsd6.prop.y+402+279*4;
P6=sum(ebsd6.isIndexed == 1);
odf6 = calcODF(ebsd6('aluminium').orientations);

grains6 = calcGrains(ebsd6('indexed'),'threshold',15*degree)
grains6 = smooth(grains6,5);
% Create border grain array
i=1;
Border6=[];
while i<=numel(grains6.boundary.grainId)/2
    if grains6.boundary.grainId(i,1)==0
        Border6(end+1)=grains6.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains6 = grains6;
selected_grains6(Border6)=[]
[omega6,a6,b6] = selected_grains6.fitEllipse;

%% SAMPLE 7

% create an EBSD variable containing the data
ebsd7 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND7.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd7.prop.y=ebsd7.prop.y+402+279*5;
P7=sum(ebsd7.isIndexed == 1);
odf7 = calcODF(ebsd7('aluminium').orientations);

grains7 = calcGrains(ebsd7('indexed'),'threshold',15*degree)
grains7 = smooth(grains7,5);
% Create border grain array
i=1;
Border7=[];
while i<=numel(grains7.boundary.grainId)/2
    if grains7.boundary.grainId(i,1)==0
        Border7(end+1)=grains7.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains7 = grains7;
selected_grains7(Border7)=[]
[omega7,a7,b7] = selected_grains7.fitEllipse;

%% SAMPLE 8

% create an EBSD variable containing the data
ebsd8 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND8.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd8.prop.y=ebsd8.prop.y+402+279*6;
P8=sum(ebsd8.isIndexed == 1);
odf8 = calcODF(ebsd8('aluminium').orientations);

grains8 = calcGrains(ebsd8('indexed'),'threshold',15*degree)
grains8 = smooth(grains8,5);
% Create border grain array
i=1;
Border8=[];
while i<=numel(grains8.boundary.grainId)/2
    if grains8.boundary.grainId(i,1)==0
        Border8(end+1)=grains8.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains8 = grains8;
selected_grains8(Border8)=[]
[omega8,a8,b8] = selected_grains8.fitEllipse;

%% SAMPLE 9

% create an EBSD variable containing the data
ebsd9 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND9.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd9.prop.y=ebsd9.prop.y+402+279*7;
P9=sum(ebsd9.isIndexed == 1);
odf9 = calcODF(ebsd9('aluminium').orientations);

grains9 = calcGrains(ebsd9('indexed'),'threshold',15*degree)
grains9 = smooth(grains9,5);
% Create border grain array
i=1;
Border9=[];
while i<=numel(grains9.boundary.grainId)/2
    if grains9.boundary.grainId(i,1)==0
        Border9(end+1)=grains9.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains9 = grains9;
selected_grains9(Border9)=[]
[omega9,a9,b9] = selected_grains9.fitEllipse;

%% SAMPLE 10

ebsd10 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_ND10.ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd10.prop.y=ebsd10.prop.y+402+279*8;
P10=sum(ebsd10.isIndexed == 1);
odf10 = calcODF(ebsd10('aluminium').orientations);

grains10 = calcGrains(ebsd10('indexed'),'threshold',15*degree)
grains10 = smooth(grains10,5);
% Create border grain array
i=1;
Border10=[];
while i<=numel(grains10.boundary.grainId)/2
    if grains10.boundary.grainId(i,1)==0
        Border10(end+1)=grains10.boundary.grainId(i,2);
    end
    
    i=i+1;
end
%Select the good grains
selected_grains10 = grains10;
selected_grains10(Border10)=[]
[omega10,a10,b10] = selected_grains10.fitEllipse;


%% FINAL ARRAYS FOR STATISTICS

TOTAL_GRAINS=numel(selected_grains1.id)+numel(selected_grains2.id)+numel(selected_grains3.id)+...
    numel(selected_grains4.id)+numel(selected_grains5.id)+numel(selected_grains6.id)+...
    numel(selected_grains7.id)+numel(selected_grains8.id)++numel(selected_grains9.id)+...
    numel(selected_grains10.id);

GRAINS.area=vertcat(selected_grains1.area,selected_grains2.area,selected_grains3.area,...
    selected_grains4.area,selected_grains5.area,selected_grains6.area,selected_grains7.area,...
    selected_grains8.area,selected_grains9.area,selected_grains10.area); %µm^2

GRAINS.area_int=round(GRAINS.area);%µm^2

GRAINS.grainSize=vertcat(selected_grains1.grainSize,selected_grains2.grainSize,...
    selected_grains3.grainSize,selected_grains4.grainSize,selected_grains5.grainSize,...
    selected_grains6.grainSize,selected_grains7.grainSize,selected_grains8.grainSize,...
    selected_grains9.grainSize,selected_grains10.grainSize); %n. of pixels

GRAINS.equivalentDiameter=2.*vertcat(selected_grains1.equivalentRadius,selected_grains2.equivalentRadius,...
    selected_grains3.equivalentRadius,selected_grains4.equivalentRadius,selected_grains5.equivalentRadius,...
    selected_grains6.equivalentRadius,selected_grains7.equivalentRadius,selected_grains8.equivalentRadius,...
    selected_grains9.equivalentRadius,selected_grains10.equivalentRadius); %µm

GRAINS.aspectRatio=1./vertcat(selected_grains1.aspectRatio,selected_grains2.aspectRatio,...
    selected_grains3.aspectRatio,selected_grains4.aspectRatio,selected_grains5.aspectRatio,...
    selected_grains6.aspectRatio,selected_grains7.aspectRatio,selected_grains8.aspectRatio,...
    selected_grains9.aspectRatio,selected_grains10.aspectRatio); % Non-dimensional

GRAINS.omega=(360/(2*pi)).*vertcat(omega1,omega2,omega3,omega4,omega5,omega6,omega7,omega8,omega9,omega10); %rad

GRAINS.meanOrientation.phi1=vertcat(selected_grains1.meanOrientation.phi1,selected_grains2.meanOrientation.phi1,...
    selected_grains3.meanOrientation.phi1,selected_grains4.meanOrientation.phi1,selected_grains5.meanOrientation.phi1,...
    selected_grains6.meanOrientation.phi1,selected_grains7.meanOrientation.phi1,selected_grains8.meanOrientation.phi1,...
    selected_grains9.meanOrientation.phi1,selected_grains10.meanOrientation.phi1); % Rad

GRAINS.meanOrientation.Phi=vertcat(selected_grains1.meanOrientation.Phi,selected_grains2.meanOrientation.Phi,...
    selected_grains3.meanOrientation.Phi,selected_grains4.meanOrientation.Phi,selected_grains5.meanOrientation.Phi,...
    selected_grains6.meanOrientation.Phi,selected_grains7.meanOrientation.Phi,selected_grains8.meanOrientation.Phi,...
    selected_grains9.meanOrientation.Phi,selected_grains10.meanOrientation.Phi); % Rad

GRAINS.meanOrientation.phi2=vertcat(selected_grains1.meanOrientation.phi2,selected_grains2.meanOrientation.phi2,...
    selected_grains3.meanOrientation.phi2,selected_grains4.meanOrientation.phi2,selected_grains5.meanOrientation.phi2,...
    selected_grains6.meanOrientation.phi2,selected_grains7.meanOrientation.phi2,selected_grains8.meanOrientation.phi2,...
    selected_grains9.meanOrientation.phi2,selected_grains10.meanOrientation.phi2); % Rad

GRAINS.TextureIndex=(textureindex(odf1)*P1+textureindex(odf2)*P2+textureindex(odf3)*P3+...
    textureindex(odf4)*P4+textureindex(odf5)*P5+textureindex(odf6)*P6+textureindex(odf7)*P7+...
    textureindex(odf8)*P8+textureindex(odf9)*P9+textureindex(odf10)*P10)/(P1+P2+P3+P4+P5+P6+P7+P8+P9+P10); 


%% ODF DATA MERGING

O1.phi1=vertcat(ebsd1('Aluminium').rotations.phi1,ebsd2('Aluminium').rotations.phi1,...
    ebsd3('Aluminium').rotations.phi1,ebsd4('Aluminium').rotations.phi1,ebsd5('Aluminium').rotations.phi1,...
    ebsd6('Aluminium').rotations.phi1,ebsd7('Aluminium').rotations.phi1,ebsd8('Aluminium').rotations.phi1,...
    ebsd9('Aluminium').rotations.phi1,ebsd10('Aluminium').rotations.phi1); % Rad

O1.Phi=vertcat(ebsd1('Aluminium').rotations.Phi,ebsd2('Aluminium').rotations.Phi,...
    ebsd3('Aluminium').rotations.Phi,ebsd4('Aluminium').rotations.Phi,ebsd5('Aluminium').rotations.Phi,...
    ebsd6('Aluminium').rotations.Phi,ebsd7('Aluminium').rotations.Phi,ebsd8('Aluminium').rotations.Phi,...
    ebsd9('Aluminium').rotations.Phi,ebsd10('Aluminium').rotations.Phi); % Rad

O1.phi2=vertcat(ebsd1('Aluminium').rotations.phi2,ebsd2('Aluminium').rotations.phi2,...
    ebsd3('Aluminium').rotations.phi2,ebsd4('Aluminium').rotations.phi2,ebsd5('Aluminium').rotations.phi2,...
    ebsd6('Aluminium').rotations.phi2,ebsd7('Aluminium').rotations.phi2,ebsd8('Aluminium').rotations.phi2,...
    ebsd9('Aluminium').rotations.phi2,ebsd10('Aluminium').rotations.phi2); % Rad

cs = crystalSymmetry('cubic');
ss = specimenSymmetry('1');

ebsd_T=orientation.byEuler(O1.phi1,O1.Phi,O1.phi2,cs,ss);
odf = calcODF(ebsd_T)
mdf = calcMDF(odf)
[density,angle] = calcAngleDistribution(mdf);
angle_deg = angle./degree;
frequency_percent = 100 * density ./sum(density);

GRAINS.TextureIndex_Good=textureindex(odf);


%%
% i=1;
% while i<=numel(GRAINS.omega)
%     if GRAINS.omega(i)>90
%         GRAINS.omega(i)=180-GRAINS.omega(i);
%     end
%     
%     i=i+1;
% end




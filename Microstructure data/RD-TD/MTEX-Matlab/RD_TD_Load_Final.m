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
setMTEXpref('xAxisDirection','north');


%% IMPORT THE DATA

%% SAMPLE 1

% create an EBSD variable containing the data
ebsd1 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_TD1 .ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
% % rotate about the z(ND)-axis
ebsd1 = rotate(ebsd1,rotation('axis',zvector,'angle',90*degree));

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
ebsd2 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_TD2 .ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd2.prop.y=ebsd2.prop.y+562;
% % rotate about the z(ND)-axis
ebsd2 = rotate(ebsd2,rotation('axis',zvector,'angle',90*degree));

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
ebsd3 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_TD3 .ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd3.prop.y=ebsd3.prop.y+562*2;
% % rotate about the z(ND)-axis
ebsd3 = rotate(ebsd3,rotation('axis',zvector,'angle',90*degree));

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
ebsd4 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_TD4 .ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd4.prop.y=ebsd4.prop.y+562*3;
% % rotate about the z(ND)-axis
ebsd4 = rotate(ebsd4,rotation('axis',zvector,'angle',90*degree));

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
ebsd5 = EBSD.load('Z:\Matlab_MTEX\mtex-5.2\RD_TD5 .ctf',CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
ebsd5.prop.y=ebsd5.prop.y+562*4;
% % rotate about the z(ND)-axis
ebsd5 = rotate(ebsd5,rotation('axis',zvector,'angle',90*degree));

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


%% FINAL ARRAYS FOR STATISTICS

TOTAL_GRAINS=numel(selected_grains1.id)+numel(selected_grains2.id)+numel(selected_grains3.id)+...
    numel(selected_grains4.id)+numel(selected_grains5.id);

GRAINS.area=vertcat(selected_grains1.area,selected_grains2.area,selected_grains3.area,...
    selected_grains4.area,selected_grains5.area); %�m^2

GRAINS.area_int=round(GRAINS.area);%�m^2

GRAINS.grainSize=vertcat(selected_grains1.grainSize,selected_grains2.grainSize,...
    selected_grains3.grainSize,selected_grains4.grainSize,selected_grains5.grainSize); %n. of pixels

GRAINS.equivalentDiameter=2.*vertcat(selected_grains1.equivalentRadius,selected_grains2.equivalentRadius,...
    selected_grains3.equivalentRadius,selected_grains4.equivalentRadius,selected_grains5.equivalentRadius); %�m

GRAINS.aspectRatio=1./vertcat(selected_grains1.aspectRatio,selected_grains2.aspectRatio,...
    selected_grains3.aspectRatio,selected_grains4.aspectRatio,selected_grains5.aspectRatio); % Non-dimensional

GRAINS.omega=(360/(2*pi)).*vertcat(omega1,omega2,omega3,omega4,omega5); %rad

GRAINS.meanOrientation.phi1=vertcat(selected_grains1.meanOrientation.phi1,selected_grains2.meanOrientation.phi1,...
    selected_grains3.meanOrientation.phi1,selected_grains4.meanOrientation.phi1,selected_grains5.meanOrientation.phi1); % Rad

GRAINS.meanOrientation.Phi=vertcat(selected_grains1.meanOrientation.Phi,selected_grains2.meanOrientation.Phi,...
    selected_grains3.meanOrientation.Phi,selected_grains4.meanOrientation.Phi,selected_grains5.meanOrientation.Phi); % Rad

GRAINS.meanOrientation.phi2=vertcat(selected_grains1.meanOrientation.phi2,selected_grains2.meanOrientation.phi2,...
    selected_grains3.meanOrientation.phi2,selected_grains4.meanOrientation.phi2,selected_grains5.meanOrientation.phi2); % Rad

GRAINS.TextureIndex=(textureindex(odf1)*P1+textureindex(odf2)*P2+textureindex(odf3)*P3+...
    textureindex(odf4)*P4+textureindex(odf5)*P5)/(P1+P2+P3+P4+P5); 

%% ODF DATA MERGING

O1.phi1=vertcat(ebsd1('Aluminium').rotations.phi1,ebsd2('Aluminium').rotations.phi1,...
    ebsd3('Aluminium').rotations.phi1,ebsd4('Aluminium').rotations.phi1,ebsd5('Aluminium').rotations.phi1); % Rad

O1.Phi=vertcat(ebsd1('Aluminium').rotations.Phi,ebsd2('Aluminium').rotations.Phi,...
    ebsd3('Aluminium').rotations.Phi,ebsd4('Aluminium').rotations.Phi,ebsd5('Aluminium').rotations.Phi); % Rad

O1.phi2=vertcat(ebsd1('Aluminium').rotations.phi2,ebsd2('Aluminium').rotations.phi2,...
    ebsd3('Aluminium').rotations.phi2,ebsd4('Aluminium').rotations.phi2,ebsd5('Aluminium').rotations.phi2); % Rad

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





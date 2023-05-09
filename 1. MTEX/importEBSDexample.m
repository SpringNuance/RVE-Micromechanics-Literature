%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('6/mmm', [5.8 5.8 4.6], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti3Al - alpha2', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('4/mmm', [4 4 4.1], 'mineral', 'TiAl - gamma', 'color', [0.56 0.74 0.56]),...
  crystalSymmetry('m-3m', [3.3 3.3 3.3], 'mineral', 'Titanium-Cubic', 'color', [0.85 0.65 0.13]),...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.94 0.5 0.5])};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files
pname = 'C:\Users\nguye\Desktop\Journal literature\Micromechanics\01_1 Exercise -MTEX';

% which files to be imported
fname = [pname '\ebsd2.ctf'];

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');


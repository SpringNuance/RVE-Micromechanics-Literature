function [F,grains] = Dream3d2DAMASKFFT(ElementNumber,MeshSize)
%% Initialize variables.
nElement_x = ElementNumber;   %%%% element number_x
nElement_y = ElementNumber;   %%%% element number_y
nElement_z = ElementNumber;   %%%% element number_z
mesh_x = MeshSize;   %%%% element number_x
mesh_y = MeshSize;   %%%% element number_y
mesh_z = MeshSize;   %%%% element number_z
filename = ('RVE_1.txt');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File name
delimiter = ' ';

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of FFT_Dream3d according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
FFT_Dream3dArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable FFT_Dream3d.
% No unimportable FFT_Dream3d rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable FFT_Dream3d, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
FFT_Dream3d = [FFT_Dream3dArray{1:end-1}];
%%%%% [Phi1] [PHI] [Phi2] [X] [Y] [Z] [Grain Id] [Phase Id]
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID FFT_Dream3dArray ans;


%% creat the grain based FFT_Dream3d
nGrain =  max(FFT_Dream3d(:,7));
F = cell(nGrain,1);
grains = zeros(nGrain,7);
grains(:,1) = 1:1:nGrain;          %%%%Grain ID
for i = 1:1:nGrain
    grains(i,6) = length(find(FFT_Dream3d(:,7)==i));
    F{i,1} = find(FFT_Dream3d(:,7)==i);
    grains(i,2) = FFT_Dream3d((F{i,1}(1,1)),8);     %%%%Grain: Phase ID
    grains(i,3) = FFT_Dream3d((F{i,1}(1,1)),1);     %%%%Grain: Euler phi1 in degree
    grains(i,4) = FFT_Dream3d((F{i,1}(1,1)),2);     %%%%Grain: Euler Phi in degree
    grains(i,5) = FFT_Dream3d((F{i,1}(1,1)),3);     %%%%Grain: Euler phi2 in degree
end
grains(:,7) = grains(:,6)/sum(grains(:,6));  %%%%Grain fraction

%% Element set for DAMASK geom file
nElement = nElement_x*nElement_y*nElement_z;
a=1:1:nElement;
b=reshape(a,nElement_x,nElement_y*nElement_z);
FFT_elementID = b';
FFT_geom = zeros(nElement_y*nElement_z,nElement_x);

for i=1:1:nElement_y*nElement_z
    for j=1:1:nElement_x
        elementID = FFT_elementID(i,j);
        FFT_geom(i,j) = FFT_Dream3d(elementID,7);
    end
end

%mesh_x_mm = mesh_x/1000000; %%%% change RVE unit: nm2mm 
%mesh_y_mm = mesh_y/1000000; %%%% change RVE unit: nm2mm 
%mesh_z_mm = mesh_z/1000000; %%%% change RVE unit: nm2mm 

%% Write DAMASK material.config - Euler angles
fn_config = sprintf('test.config');
fidw_config = fopen(fn_config,'w+');
%print materials setting
fprintf(fidw_config,'#-------------------#\r\n');    
fprintf(fidw_config,'<microstructure>\r\n');   
fprintf(fidw_config,'#-------------------#\r\n');   
for i=1:1:nGrain
    fprintf(fidw_config,'[Al %i]\r\n',i);
    fprintf(fidw_config,'crystallite 1 \r\n');  
    fprintf(fidw_config,'(constituent) phase %i texture %i fraction 1.0\r\n\r\n',grains(i,2),i);
end
fprintf(fidw_config,'#-------------------#\r\n');    
fprintf(fidw_config,'<texture>\r\n');   
fprintf(fidw_config,'#-------------------#\r\n');  
for i=1:1:nGrain
    fprintf(fidw_config,'[%i]\r\n',i); 
    fprintf(fidw_config,'(gauss)	phi1	%5.2f	Phi	%5.2f	phi2	%5.2f	scatter	0.000 fraction	1.000\r\n\r\n',grains(i,3:5));
end
fclose(fidw_config);

%% Write DAMASK geom
fn_geom = sprintf('test.geom');
fidw_geom = fopen(fn_geom,'w+');
%print materials setting
fprintf(fidw_geom,'5 header\r\n');    
fprintf(fidw_geom,'grid	a %3.i	b %3.i	c %3.i\r\n',nElement_x,nElement_y,nElement_z);  
fprintf(fidw_geom,'size	x %7.i	y %7.i	z %7.i\r\n', mesh_x, mesh_y, mesh_z);
fprintf(fidw_geom,'origin	x 0	    y 0	    z 0\r\n');  
fprintf(fidw_geom,'microstructures %i\r\n', nGrain); 
fprintf(fidw_geom,'homogenization	1\r\n');  
formatGeom = repmat('%6.i',1,nElement_x); 
for i=1:1:nElement_y*nElement_z
    fprintf(fidw_geom,formatGeom,FFT_geom(i,:));
    fprintf(fidw_geom,'\r\n');  
end
fclose(fidw_geom);
end


% Points for mesh generation in GMSH

clear all;
close all;


% Example with cell
I = imread('Images/flip_05_singlecell.png');
phi = chanvese(I,1000,350,1,40,1000);

[C,h] = contour(phi, [0 0], 'r', 'LineWidth',2);    %make points
C = C(1:2,2:(length(C)-1));                         %Points to vector
C(2,1:length(C)) = size(phi,1)-C(2,1:length(C));    %Turn y-axis

k = 5;
C_new = C(1:2,1:k:(length(C)));                         %Take only every k'th point
contour(phi, [0 0], 'r', 'LineWidth',2);            %plot
scatter(C_new(1,1:length(C_new)),C_new(2,1:length(C_new)));         %plot




n = 1:length(C_new);
c = [n;C_new];
m = [n(1:length(n)-1);n(1:length(n)-1);n(2:length(n))];


% Output for GMSH
fileID = fopen('contourpoints.geo','w');
fprintf(fileID, 'lc = 1e+500; \n'); 
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc }; \n', c);
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', m);
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', [length(n);length(n);1]);
fprintf(fileID,'Line Loop(%d) = {', length(n)+1);
fprintf(fileID,'%d,', 1:(length(C_new)-1));
fprintf(fileID,'%d}; \n', length(C_new));
fprintf(fileID,'Plane Surface(%d) = {%d} ;', [length(C_new)+2,length(C_new)+1]);
fclose(fileID);





% Example for nucleus
I = imread('Images/flip_05_nucleus.png');
phi = chanvese(I,1000,350,1,1,1000);

[C,h] = contour(phi, [0 0], 'r', 'LineWidth',2);    %make points
C = C(1:2,2:(length(C)-1));                         %Points to vector
C(2,1:length(C)) = size(phi,1)-C(2,1:length(C));    %Turn y-axis

k = 3;
C = C(1:2,1:k:(length(C)));                         %Take only every k'th point

contour(phi, [0 0], 'r', 'LineWidth',2);            %plot
scatter(C(1,1:length(C)),C(2,1:length(C)));         %plot


% Coordinate output
n = 1:length(C);
c = [n;C];
fileID = fopen('subdomain.txt','w');
fprintf(fileID,'%f,%f\n', c(2:3,n));
fclose(fileID);








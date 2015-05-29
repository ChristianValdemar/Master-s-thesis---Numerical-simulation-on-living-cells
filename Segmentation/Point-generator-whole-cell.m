% Points for mesh generation in GMSH
% This generator generates a membrane and scales the output points, such
% that they are in the propper scale. 
% It adds "dummy" points near the membrane to avoid to fine mesh in the rest of the cell. 

clear all;
close all;

scalar = 0.0448; %um/pixel
mem_thickness = 0.009; %um
mem = mem_thickness/scalar; % Membrane tickness in pixel





% First make cell boundary
I = imread('Images/flip_05_singlecell.png');
phi_org = chanvese(I,1000,1000,1,40,1000);                   %Levelset func.
phi = flipud(phi_org);                                      %Flip y-axis 

[C,h] = contour(phi, [0 0], 'r', 'LineWidth',2);            %make points
C = C(1:2,2:(length(C)-1));                                 %Points to vector


k = 6;
C_cell = C(1:2,1:k:(length(C)));                            %Take only every k'th point
contour(phi, [0 0], 'r', 'LineWidth',2);                    %plot
scatter(C_cell(1,1:length(C_cell)),C_cell(2,1:length(C_cell))); %plot










%
% Make boundary for nucleus
%

I = imread('Images/flip_05_nucleus.png');
phi_nuc_org = chanvese(I,1000,350,1,5,1000);
phi_nuc = flipud(phi_nuc_org);                              %Flip y-axis 
phi_norm = -Normal_phi(phi_nuc);                             %Unit normals for phi


[C_n,h] = contour(phi_nuc, [0 0], 'r', 'LineWidth',2);      %make points
C_n = C_n(1:2,2:(length(C_n)-1));                           %Points to vector

k = 2;
C_nuc = C_n(1:2,1:k:(length(C_n)));                         %Take only every k'th point

contour(phi, [0 0], 'r', 'LineWidth',2);                    %plot
scatter(C_nuc(1,1:length(C_nuc)),C_nuc(2,1:length(C_nuc))); %plot
%hold on




%
% Make boundary bleaching area
%

I = imread('Images/flip_05_eGFP-withROI2.tif');
phi_bleach_org = chanvese(I,1000,1,1,100,1000,1);
phi_bleach = flipud(phi_bleach_org);                           %Flip y-axis 


[C_b,h] = contour(phi_bleach, [0 0], 'r', 'LineWidth',2);   %make points
C_b = C_b(1:2,2:(length(C_b)-1));                           %Points to vector

k = 2;
C_bl = C_b(1:2,1:k:(length(C_b)));                          %Take only every k'th point

contour(phi_bleach, [0 0], 'r', 'LineWidth',2);             %plot
scatter(C_bl(1,1:length(C_bl)),C_bl(2,1:length(C_bl)));     %plot
%hold on



scatter(C_cell(1,1:length(C_cell)),C_cell(2,1:length(C_cell))); %plot
hold on
scatter(C_nuc(1,1:length(C_nuc)),C_nuc(2,1:length(C_nuc))); %plot
scatter(C_bl(1,1:length(C_bl)),C_bl(2,1:length(C_bl)));     %plot


%
% Clear variables
%

clearvars -except C_nuc C_cell C_n phi_nuc_org phi_nuc phi_norm phi_org phi C scalar mem mem_thickness C_bl C_b phi_bleach_org phi_bleach saved_cell saved_nuc saved_bleach


%
% Take the full cell and nucleus
%
center = [(max(C_nuc(1,:))+min(C_nuc(1,:)))/2;(max(C_nuc(2,:))+min(C_nuc(2,:)))/2];
P_nuc = C_nuc;
P_cell = C_cell;


%
% Make membrane
%

% Move points in outward normal direction for membrane
x = P_nuc(1,:);
y = P_nuc(2,:);
index1 = sub2ind(size(phi_norm),round(y),round(x),ones(size(x)));
index2 = sub2ind(size(phi_norm),round(y),round(x),ones(size(x))*2);
x1 = x + 1/2*mem*phi_norm(index1);
y1 = y + 1/2*mem*phi_norm(index2);
C_big(1,:) = x1;
C_big(2,:) = y1;
%plot(C_big(1,:),C_big(2,:),'.-');


% Move points in inward normal direction for membrane
x1 = x - 1/2*mem*phi_norm(index1);
y1 = y - 1/2*mem*phi_norm(index2);
C_small(1,:) = x1;
C_small(2,:) = y1;
%plot(C_small(1,:),C_small(2,:),'r.-');


% Move points in inward normal direction for mesh refinement
k = 3;
C_nuc1 = P_nuc(1:2,2:k:(length(P_nuc))-1);                         %Take only every k'th point
x = C_nuc1(1,:);
y = C_nuc1(2,:);
index1 = sub2ind(size(phi_norm),round(y),round(x),ones(size(x)));
index2 = sub2ind(size(phi_norm),round(y),round(x),ones(size(x))*2);
x1 = x - 10*mem*phi_norm(index1);
y1 = y - 10*mem*phi_norm(index2);
C_small2(1,:) = x1;
C_small2(2,:) = y1;

% Move points in outward normal direction for mesh refinement
x1 = x + 10*mem*phi_norm(index1);
y1 = y + 10*mem*phi_norm(index2);
C_big2(1,:) = x1;
C_big2(2,:) = y1;


% Move points in outward normal direction for initilization of membrane
% Move the points a bit more than the thickness of the membrane 
x = P_nuc(1,:);
y = P_nuc(2,:);
index1 = sub2ind(size(phi_norm),round(y),round(x),ones(size(x)));
index2 = sub2ind(size(phi_norm),round(y),round(x),ones(size(x))*2);
x1 = x + 0.6*mem*phi_norm(index1);
y1 = y + 0.6*mem*phi_norm(index2);
C_init_out(1,:) = x1;
C_init_out(2,:) = y1;

% Move points in inward normal direction for initilization of membrane
% Move the points a bit more than the thickness of the membrane 
x1 = x - 0.6*mem*phi_norm(index1);
y1 = y - 0.6*mem*phi_norm(index2);
C_init_in(1,:) = x1;
C_init_in(2,:) = y1;


%
% Scale data 
%

S_cell  = P_cell.*scalar;
S_nuc   = P_nuc.*scalar;
S_small = C_small.*scalar;
S_small2= C_small2.*scalar;
S_big   = C_big.*scalar;
S_big2  = C_big2.*scalar;
center  = center.*scalar;
S_bl  = C_bl.*scalar;
S_init_out  = C_init_out.*scalar;
S_init_in  = C_init_in.*scalar;


%
% Output to file
%


% Coordinates for nucleus
n_nuc = 1:length(S_nuc);
c_nuc = [n_nuc;S_nuc];
fileID = fopen('Mesh/nucleus.txt','w');
fprintf(fileID,'%f,%f\n', c_nuc(2:3,n_nuc));
fclose(fileID);

% Coordinates for the inner wall of the membrane
n_small = 1:length(S_small);
c_small = [n_small;S_small];
fileID = fopen('Mesh/membrane_in.txt','w');
fprintf(fileID,'%f,%f\n', c_small(2:3,n_small));
fclose(fileID);

% Coordinates for the outer wall of the membrane
n_big = 1:length(S_big);
c_big = [n_big;S_big];
fileID = fopen('Mesh/membrane_out.txt','w');
fprintf(fileID,'%f,%f\n', c_big(2:3,n_big));
fclose(fileID);

% Coordinates for bleaching area
n_bl = 1:length(S_bl);
c_bl = [n_bl;S_bl];
fileID = fopen('Mesh/bleaching_area.txt','w');
fprintf(fileID,'%f,%f\n', c_bl(2:3,n_bl));
fclose(fileID);

% Coordinates for the outer wall of the membrane for initilization
n_init_out = 1:length(S_init_out);
c_init_out = [n_init_out;S_init_out];
fileID = fopen('Mesh/init_membrane_out.txt','w');
fprintf(fileID,'%f,%f\n', c_init_out(2:3,n_init_out));
fclose(fileID);

% Coordinates for the inner wall of the membrane for initilization
n_init_in = 1:length(S_init_in);
c_init_in = [n_init_in;S_init_in];
fileID = fopen('Mesh/init_membrane_in.txt','w');
fprintf(fileID,'%f,%f\n', c_init_in(2:3,n_init_in));
fclose(fileID);



% Without nucleus edge, only membrane walls

% Output for GMSH
% Cell boundary, membrane and nucleus boundary
n = 1:length(S_cell);
c_cell = [n;S_cell];
m = [n(1:length(n)-1);n(1:length(n)-1);n(2:length(n))];

n_start1 = length(S_cell)+1;
n_end1 = length(S_cell)+length(S_small);
n_small = (n_start1:n_end1);
c_small = [n_small;S_small];
m_small = [n_small(1:length(n_small)-1);n_small(1:length(n_small)-1);n_small(2:length(n_small))];

n_start2 = n_start1 + length(S_small);
n_end2 = n_end1+length(S_big);
n_big = (n_start2:n_end2);
c_big = [n_big;S_big];
m_big = [n_big(1:length(n_big)-1);n_big(1:length(n_big)-1);n_big(2:length(n_big))];

n_start3 = n_start2 + length(S_big);
n_end3 = n_end2+length(S_small2);
n_small2 = (n_start3:n_end3);
c_small2 = [n_small2;S_small2];
m_small2 = [n_small2(1:length(n_small2)-1);n_small2(1:length(n_small2)-1);n_small2(2:length(n_small2))];

n_start4 = n_start3 + length(S_small2);
n_end4 = n_end3+length(S_big2);
n_big2 = (n_start4:n_end4);
c_big2 = [n_big2;S_big2];
m_big2 = [n_big2(1:length(n_big2)-1);n_big2(1:length(n_big2)-1);n_big2(2:length(n_big2))];

n_start5 = n_start4 + length(S_big2);
n_end5 = n_end4+length(S_bl);
n_bl = (n_start5:n_end5);
c_bl = [n_bl;S_bl];
m_bl = [n_bl(1:length(n_bl)-1);n_bl(1:length(n_bl)-1);n_bl(2:length(n_bl))];

n_start6 = n_start5 + length(S_bl);
n_end6 = n_end5+length(S_nuc);
n_nuc = (n_start6:n_end6);
c_nuc = [n_nuc;S_nuc];
m_nuc = [n_nuc(1:length(n_nuc)-1);n_nuc(1:length(n_nuc)-1);n_nuc(2:length(n_nuc))];


%Print with extra line in membrane
fileID = fopen('Mesh/contourpoints.geo','w');
fprintf(fileID, 'lc1 = 0.5; \n');           % Mesh size for cell boundary
fprintf(fileID, 'lc2 = 0.002; \n');          % Mesh size for inner wall in membrane
fprintf(fileID, 'lc3 = 0.002; \n');          % Mesh size for outer wall in membrane
fprintf(fileID, 'lc5 = 0.3; \n');           % Mesh size for inner mesh refinement
fprintf(fileID, 'lc6 = 0.3; \n');           % Mesh size for outer mesh refinement
fprintf(fileID, 'lc7 = 0.3; \n');           % Mesh size for center point in nucleus
fprintf(fileID, 'lc8 = 0.002; \n');          % Mesh size for line in membrane
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc1 }; \n', c_cell);     % Cell boundary points
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc2 }; \n', c_small);    % Membrane - Inner wall
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc3 }; \n', c_big);      % Membrane - Outer wall
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc5 }; \n', c_small2);   % Inner refinement points
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc6 }; \n', c_big2);     % Outer refinement points
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc1 }; \n', c_bl);     % Bleaching area
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc8 }; \n', c_nuc);     % Extra line in membrane
fprintf(fileID,'Point(%d) = { %f , %f ,0 , lc7 };  \n', [n_end6+1;center]);      % Center point for mesh refinement
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', m);                       % Lines for cell bound
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', [length(n);length(n);1]); 
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', m_small);                 % Lines for Inner wall
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', [n_end1;n_end1;n_start1]);
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', m_big);                   % Lines for Outer wall
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', [n_end2;n_end2;n_start2]);
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', m_bl);                    % Bleaching area
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', [n_end5;n_end5;n_start5]);
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', m_nuc);                   % Extra line in nucleus
fprintf(fileID,'Line(%d) = {%d,%d} ; \n', [n_end6;n_end6;n_start6]);
fprintf(fileID,'Line Loop(%d) = {', n_end3+1);                      % Line loop for cell bound
fprintf(fileID,'%d,', 1:(length(S_cell)-1));   
fprintf(fileID,'%d}; \n', length(S_cell));
fprintf(fileID,'Line Loop(%d) = {', n_end3+2);                      % Line loop for Inner wall
fprintf(fileID,'%d,', n_start1:(n_end1-1));
fprintf(fileID,'%d}; \n', n_end1);
fprintf(fileID,'Line Loop(%d) = {', n_end3+3);                      % Line loop for Outer wall
fprintf(fileID,'%d,', n_start2:(n_end2-1));
fprintf(fileID,'%d}; \n', n_end2);
fprintf(fileID,'Line Loop(%d) = {', n_end3+4);                      % Line loop for Bleaching area
fprintf(fileID,'%d,', n_start5:(n_end5-1));
fprintf(fileID,'%d}; \n',n_end5);
fprintf(fileID,'Line Loop(%d) = {', n_end3+5);                      % Line loop for extra line in membrane
fprintf(fileID,'%d,', n_start6:(n_end6-1));
fprintf(fileID,'%d}; \n',n_end6);
fprintf(fileID,'Plane Surface(%d) = {%d,%d,%d}; \n', [3,n_end3+1,n_end3+3,n_end3+4]); % Cytoplasma
fprintf(fileID,'Plane Surface(%d) = {%d,%d}; \n', [0,n_end3+3,n_end3+5]); % Outer Membrane
fprintf(fileID,'Plane Surface(%d) = {%d,%d}; \n', [1,n_end3+5,n_end3+2]); % Inner Membrane
fprintf(fileID,'Plane Surface(%d) = {%d}; \n', [4,n_end3+2]); % Nucleus
fprintf(fileID,'Plane Surface(%d) = {%d}; \n', [2,n_end3+4]); % Bleaching area
fprintf(fileID,'Point{%d} In Surface{%d}; \n', [n_end6+1,4]);
fprintf(fileID,'Point{%d} In Surface{%d}; \n', [c_small2(1,:);(4)*ones(size(c_small2(1,:)))]);        % Mesh refinement points in nucleus
fprintf(fileID,'Point{%d} In Surface{%d}; \n', [c_big2(1,:);3*ones(size(c_big2(1,:)))]);          % Mesh refinement points in cytop.
fclose(fileID);

























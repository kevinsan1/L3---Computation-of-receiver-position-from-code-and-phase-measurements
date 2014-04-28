%% Lab 3
% Difference Equations for GPS
clear all; clc;
%% Import file
addpath('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L3 - Computation of receiver position from code and phase measurements')
excfile = xlsread('DiffAssignmentL1');
excfile = padarray(excfile,5,'pre');
cellDiff = num2cell(excfile); % same as a cell
%% Constants
c = 299792458; % speed of light (m/s)
[lambda1,lambda2] = cellDiff{78:79,2};
mu = 3.986005e14; % universal gravitational parameter (m/s)^3
omega_e_dot = 7.2921151467e-5; % earth rotation rate (rad/s)
F = -4.442807633e-10; % s/m^1/2
nRows = 5;
%% Variables
dtA = excfile(24,3);
dtB = excfile(25,3);
sNum = excfile((28:33),1);
refSatNum = 20;
for i = 1:length(sNum) % vector of ones at satellite numbers
    vectorOfOnes(sNum(i),1) = 1; 
end
[Xref,Xrov,Yref,Yrov,Zref,Zrov] = cellDiff{(37:38),3:5};
for i = 43:48 % 
    rhoDot(excfile(i-15,1),1) = excfile(i-15,2);
    satCoordinates(excfile(i,1),:) = excfile(i,3:5);
    rhoAtoP(excfile(i+23,1),:) = excfile(i+23,2);
    rhoB0toP(excfile(i+23,1),:) = excfile(i+23,3);
end
count = 1;
for i = 6:3:21
   P1ref(sNum(count),1) = excfile(i,3);
   P1rov(sNum(count),1) = excfile(i+1,3);
   phi1ref(sNum(count),1) = excfile(i,4);
   phi1rov(sNum(count),1) = excfile(i+1,4);
   count = count + 1;
end
weightMatrix = excfile(52:61,1:10);
%% Main steps of Calculation
%% 1. Synchronize observables from both receivers using equations
% 14
P1refeq14 = P1ref + rhoDot*dtA;
P1roveq14 = P1rov + rhoDot*dtB;
phi1refeq14 = lambda1*phi1ref + rhoDot * dtA;  
phi1roveq14 = lambda1*phi1rov + rhoDot * dtB;  
%% 2. Compute single and double differences - equations (15) and
% (18). Use satellite 20 as a reference satellite
% for double differencing.
P1eq15 = P1refeq14 - P1roveq14;
phi1eq15 = (phi1refeq14 - phi1roveq14);
dtAB = dtB - dtA;
% Double differences
phi1eq18 = phi1eq15 - phi1eq15(refSatNum)*vectorOfOnes;
P1eq18 = P1eq15 - P1eq15(refSatNum)*vectorOfOnes;
%% 3. Compute coefficients a , a , a and ?pq - equations (12).
% XYZ
% ref
ax_S_rov = -1 * (satCoordinates(:,1)-Xrov*vectorOfOnes)./(rhoB0toP);
ay_S_rov = -1 * (satCoordinates(:,2)-Yrov*vectorOfOnes)./(rhoB0toP);
az_S_rov = -1 * (satCoordinates(:,3)-Zrov*vectorOfOnes)./(rhoB0toP);
% rov
% ax_S_rov = -1 * (satCoordinates(:,1)-Xrov*vectorOfOnes)./(rhoB0toP);
% ay_S_rov = -1 * (satCoordinates(:,2)-Yrov*vectorOfOnes)./(rhoB0toP);
% az_S_rov = -1 * (satCoordinates(:,3)-Zrov*vectorOfOnes)./(rhoB0toP);
%
aXB_st = -ax_S_rov + ax_S_rov(20);
aYB_st = -ay_S_rov + ay_S_rov(20);
aZB_st = -az_S_rov + az_S_rov(20);
rhoAB_0 = rhoAtoP - rhoB0toP;
%% 4. Fill in matrixes A, L and compute least square solution of
% equations (21).
p_AB_s = rhoB0toP-rhoAtoP;
nv = [4,5,6,24,25];
A = [aXB_st(nv,1), aYB_st(nv,1), aZB_st(nv,1)]; 
A = [ A , eye(length(nv)) * lambda1 ; A , zeros(length(nv)) ];
% L
L = [phi1eq18(nv) - rhoAB_0(nv) + ones(5,1)*rhoAB_0(20); P1eq18(nv) - rhoAB_0(nv) + ones(5,1)*rhoAB_0(20)];
%
X = inv(A'*weightMatrix*A)*A'*weightMatrix*L;

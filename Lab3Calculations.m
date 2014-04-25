%% Lab 3
% Difference Equations for GPS
clear all; clc;
%% Import file
addpath('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L3 - Computation of receiver position from code and phase measurements')
mmainFile = xlsread('DiffAssignmentL1');
mmainFile = padarray(mmainFile,5,'pre');
cellDiff = num2cell(mmainFile); % same as a cell
%% Constants
c = 299792458; % speed of light (m/s)
[lambda1,lambda2] = cellDiff{78:79,2};
mu = 3.986005e14; % universal gravitational parameter (m/s)^3
omega_e_dot = 7.2921151467e-5; % earth rotation rate (rad/s)
F = -4.442807633e-10; % s/m^1/2
nRows = 5;
%% Variables
dtA = mmainFile(24,3);
dtB = mmainFile(25,3);
satelliteNumbers = mmainFile((28:33),1);
[Xref,Xrov,Yref,Yrov,Zref,Zrov] = cellDiff{(37:38),3:5};
for i = 43:48 % 
    rhoDot(mmainFile(i-15,1),1) = mmainFile(i-15,2);
    satCoordinates(mmainFile(i,1),:) = mmainFile(i,3:5);
    rhoAtoP(mmainFile(i+23,1),:) = mmainFile(i+23,2);
    rhoB0toP(mmainFile(i+23,1),:) = mmainFile(i+23,3);
end
count = 1;
for i = 6:3:21
   P1ref(satelliteNumbers(count),1) = mmainFile(i,3);
   P1rov(satelliteNumbers(count),1) = mmainFile(i+1,3);
   phi1ref(satelliteNumbers(count),1) = mmainFile(i,4);
   phi1rov(satelliteNumbers(count),1) = mmainFile(i+1,4);
   count = count + 1;
end
weightMatrix = mmainFile(52:61,1:10);
%% Main steps of Calculation
%% 1. Synchronize observables from both receivers using equations
% 14
P1refeq14 = P1ref + rhoDot*dtA;
P1roveq14 = P1rov + rhoDot*dtB;
phi1refeq14 = lambda1*phi1ref + rhoAtoP * dtA; % lambda times phi
phi1roveq14 = lambda1*phi1rov + rhoAtoP * dtB; % lamda times phi
%% 2. Compute single and double differences ? equations (15) and
% (18). Use satellite 20 as a
% AB,0
% i dd
% reference satellite for double differencing.

%% 3. Compute coefficients a , a , a and ?pq - equations (12).
% XYZ

%% 4. Fill in matrixes A, L and compute least square solution of
% equations (21).

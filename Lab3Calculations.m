%% Lab 3
% Difference Equations for GPS
clear all; clc;
%% Import file
addpath('/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L3 - Computation of receiver position from code and phase measurements')
DiffAssignmentL1 = xlsread('DiffAssignmentL1');
cellDiff = num2cell(DiffAssignmentL1); % same as a cell
%% Constants
c = 299792458; % speed of light (m/s)
mu = 3.986005e14; % universal gravitational parameter (m/s)^3
omega_e_dot = 7.2921151467e-5; % earth rotation rate (rad/s)
F = -4.442807633e-10; % s/m^1/2
nLineUpRows = 5;
%% Variables
P1 = DiffAssignmentL1((6:22)-nLineUpRows,3);
phi1 = DiffAssignmentL1((6:22)-nLineUpRows,4);
dtA = DiffAssignmentL1(24-nLineUpRows,3);
dtB = DiffAssignmentL1(25-nLineUpRows,3);
rhoDot = DiffAssignmentL1((28:33)-nLineUpRows,2);
satelliteNumbers = DiffAssignmentL1((28:33)-nLineUpRows,1);
[Xref,Xrov,Yref,Yrov,Zref,Zrov] = cellDiff{(37:38)-nLineUpRows,3:5};
%% Satellite Coordinates
sat(20,:) = DiffAssignmentL1(43-nLineUpRows,3:5);
sat(4,:) = DiffAssignmentL1(44-nLineUpRows,3:5);
sat(5,:) = DiffAssignmentL1(45-nLineUpRows,3:5);
sat(6,:) = DiffAssignmentL1(46-nLineUpRows,3:5);
sat(24,:) = DiffAssignmentL1(47-nLineUpRows,3:5);
sat(25,:) = DiffAssignmentL1(48-nLineUpRows,3:5);


%% Main steps of Calculation
%% 1. Synchronize observables from both receivers using equations
% (1
%% 2. Compute single and double differences ? equations (15) and
% (18). Use satellite 20 as a
% AB,0
% i dd
% reference satellite for double differencing.

%% 3. Compute coefficients a , a , a and ?pq - equations (12).
% XYZ

%% 4. Fill in matrixes A, L and compute least square solution of
% equations (21).

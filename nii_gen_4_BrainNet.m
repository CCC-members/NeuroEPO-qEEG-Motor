% Script to generate .nii file to be loaded in BrainNet Viewer
% for visualization of CCA scores on Cortex
% Authors: Usama Riaz, Fuleah A. Razzaq, Pedro A Valdes Sosa
clear Y th1;
clc;
load('CogEEG.mat');% load cognition CCA scores
CogEEG = reshape(CogEEG,[3244 49]);%reshape into a 2-D matrix where each column represents a frequency
freq = 0.3906:0.3906:19.15;% 49 frequencies
delta = 4:9;%defining delta band
theta = 9:19;%defining theta band
alpha = 19:32;%defining alpha band
beta = 34:49;%defining beta band
Y = log((10000*((max(MotorEEG(:,theta),[],2)))).^2);% data transformation for visualizatio clarity
th1 = prctile(Y,95);%defining threshold for viewing onl 5% elements
Y(Y<th1)=-120;% appying threshol
max(Y(Y>th1))%maximum colormap value for a specific band
min(Y(Y>th1))%minimum colormap value for a specific band
XYZ_to_AAL_BrainNet_coord(Y, "MotorCCA_delta.nii");% saving te .nii file for BrainNet Viewer

clc
clear
currentFolder = pwd;
folders = genpath(currentFolder);
addpath(folders);
load("iris_data.mat")

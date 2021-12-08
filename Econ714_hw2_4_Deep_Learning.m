% -------------------------------------------------------------------------
% Econ 714
% Homework 2
% Fall 2021 
% Javier Tasso
% Deep Learning
% -------------------------------------------------------------------------
clearvars
clc 
cd 'C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2'
% -------------------------------------------------------------------------

% Define parameters 
gamma = 2; 
rho = 0.04;
A = 0.5;
alpha = 0.36;
delta = 0.05;
k_min = 0.1;
k_max = 10;
batch_size = 1000;
n_burned = 1;
n_neurons = 8;
input_dim = 1;
dt = 0.1;
epsilon = 1; 

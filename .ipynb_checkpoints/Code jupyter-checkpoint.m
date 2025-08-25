% ---
% jupyter:
%   jupytext:
%     formats: ipynb,m:percent
%     text_representation:
%       extension: .m
%       format_name: percent
%       format_version: '1.3'
%       jupytext_version: 1.17.2
%   kernelspec:
%     display_name: Matlab
%     language: matlab
%     name: matlab
% ---

% %% [markdown]
% # A Finite Element Solution for a plate with MATLAB code.

% %% [markdown]
% The problem is show in the figure 6-16. 

% %% [markdown]
% ![Figure 1](figure1.jpg)

% %% [markdown]
% This problem is show in the chapter 6.5 of the book: *A firts Course in the Finite Element Method* written by Daryl L. Logan. In this 
% chapter a solution is make by subdivided the plate in two CST. The purpose of this project is to make more divisions considering 400 CST and see the results. 

% %% [markdown]
% ## Preprocesing 

% %%
% Geometry definition

% The values are: 

t = 1; % in
E = 30*10^6; %psi
v = 0.3; 

%The mesh is divided in 200 squares, 10 division in the y direction and 20 divisions in the x direction. With this division each rectangle size 1x1 in^2. 
%At each rectangle the CST is generated joining the lower-left node with the upper-right node of the rectangle. 
%The vectors for this divisions are defined: 

x0 = 0;
y0 = 0;
x1 = 20;
y1 = 10; 
x = linspace(x0,x1,20); 
y = linspace(y0,y1,10); 
division = zeros(length(x),length(y));


% %%

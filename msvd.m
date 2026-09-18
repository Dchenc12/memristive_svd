clear
close all
run("init.m")


%%
clc

voltage = 0.2; % V
Gmin = 0;   % uS
Gmax = 30;  % uS
Wr_E = 0.35; % uS
miu_Rd_E = 0.0;% uS
sig_Rd_E = 0.3;% uS

vrpt = 4;
block_size = 0;

sizeX = size(dataX);
singular_num = sizeX(2);
dataX_max = max(dataX,[],"all");


Memris = [];
division = [2,5];
outer_iter = 10;

seed = 2;


[U,S,V] = svd(dataX);
S = diag(S);
sign_first = sign(V(1,:));
sign_first_U = sign(U(1,:));

rng(seed);

fprintf('V ground truth is: \n')
disp(V)
fprintf('S ground truth is: \n')
disp(S')
fprintf('Solving with SREA\n')
[V_w_SREA,S_w_SREA,U_w_SREA] = svd_decomposition_w_srea_(dataX,division,Gmin,Gmax,Wr_E,miu_Rd_E,sig_Rd_E,voltage,block_size,vrpt,20);

V_eigen_q_n_orig = V_w_SREA;
V_w_SREA = get_direction(V,V_w_SREA);

fprintf('V solved with SREA is: \n')
disp(V_w_SREA)
fprintf('S solved with SREA is: \n')
disp(S_w_SREA)

division_analog = 16;
vrpt_wo_srea = 1;
fprintf('Solving without SREA\n')

[V_wo_SREA,S_wo_SREA,U_wo_SREA] = svd_decomposition_wo_srea_(dataX,division_analog,Gmin,Gmax,Wr_E,miu_Rd_E,sig_Rd_E,voltage,block_size,vrpt_wo_srea,20);
V_eigen_n_nspa_orig = V_wo_SREA;
V_wo_SREA = get_direction(V,V_wo_SREA);
fprintf('V solved without SREA is: \n')
disp(V_wo_SREA)
fprintf('S solved without SREA is: \n')
disp(S_wo_SREA)
        





















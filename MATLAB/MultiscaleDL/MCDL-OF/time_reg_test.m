[y,y_true,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(0,'dissertation_adjust2');
K = 2;
J = 8;

kernel = [[ 0  0  0;
            0  1  0;
            0  0  0],...
          [-1 -1 -1;
           -1 -1 -1;
           -1 -1 -1]];

compute_time_reg(kernel,Xtrue,K,J)

kernel = [[ 1  1  1;
            1  1  1;
            1  1  1],...
          [-1 -1 -1;
           -1 -1 -1;
           -1 -1 -1]];

compute_time_reg(kernel,Xtrue,K,J)

kernel = [[ 1],...
          [-1]];

compute_time_reg(kernel,Xtrue,K,J)

[R, dR_flat] = compute_time_reg_softmin_flat(Xtrue, K, J, 1e-3);
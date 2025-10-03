function [Xf, cgst, opt] = x_update_softmin_lbfgs_gpu(X0,Df,bf,Sf,Y,U,opt,N,M,K,J,T)
    if opt.lambda2 > 0
        if opt.ism_init
            Xf0 = solvedbi_sm(Df, opt.rho1, bf);
            X0 = ifft2(Xf0,'symmetric');
        end
        
        % LBFGS Optimizer
        options = optimoptions('fminunc', ...
            'Algorithm','quasi-newton', ...        % this is L-BFGS in modern MATLAB
            'HessUpdate','lbfgs', ...              % force limited-memory BFGS
            'SpecifyObjectiveGradient',true, ...
            'MaxIterations',200, ...
            'Display','iter', ...
            'OptimalityTolerance',1e-10, ...
            'StepTolerance',1e-12);

        [Xvec,~,~,output] = fminunc(@(x)obj_and_grad_x(x,Df,Sf,Y,U,opt.rho1,opt.lambda2,K,J), X0(:), options);
        cgst.pit = output.iterations;
        X = reshape(Xvec,size(X0));
        Xf = fft2(X);
    else
        Xf = solvedbi_sm(Df, opt.rho1, bf);
        cgst.pit = 1;
    end

end
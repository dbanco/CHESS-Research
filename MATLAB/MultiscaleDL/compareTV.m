% Compare the TV terms
sets = [13,15,14];
fprintf('\lambda_2     TV    ||x||_1    ')
for i = 1:numel(sets)
    load(['C:\Users\dpqb1\Documents\Outputs\',...
         sprintf('toy1_exp_TV%i_sig_5',sets(i)),'\output_18.mat']);
    
    DPY = DiffPhiX_1D(outputs.X);  
    Jtv = sum(abs(vec(DPY)));
    Jl1 = sum(abs(vec(bsxfun(@times, opt.L1Weight, Y))));
    pVars = [outputs.lambda1 outputs.lambda2 Jtv Jl1];
    sprintf('%0.3f     %0.3f     %0.3f     %0.3f\n',pVars)

end

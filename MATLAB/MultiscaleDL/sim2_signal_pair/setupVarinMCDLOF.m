function varin = setupVarinMCDLOF(opt,topDir,dataset,K,scales,sig_ind,lambdaVals,lambdaOFVals,lambdaHSVals,ind1,ind2,ind3)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for sig_i = sig_ind
    % j_s_select = find(lambdaVals == selected_lam_s_vec(sig_i));
    for j_s = ind1
        for j_hs = ind2
            for j_of = ind3
                varin = {lambdaVals,lambdaOFVals,lambdaHSVals,...
                        j_s,j_of,j_hs,sigmas,sig_i,opt,topDir,dataset,K,scales};
                save(fullfile(jobDir,['varin_',num2str(k),'.mat']),'varin','funcName')
                k = k + 1;
            end
        end
    end
end

end
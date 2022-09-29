function figOut = waterfallFigure(nnInd,MM,top_dir,sim_name,ps_name,P)
alg_name = 'CoupledCGTV';
coupledPS_dir = fullfile(top_dir,sim_name,alg_name);

figOut = figure;
i = 1;
for nn = nnInd
    % Load param search data
    try
        c_data = load(fullfile(coupledPS_dir,[ps_name,'_',num2str(nn)]));
        [N,T] = size(c_data.B);
        A0 = unshifted_basis_vector_ft_stack_zpad(c_data.P);
        fits = Ax_ft_1D_Time(A0,...
            squeeze(c_data.X_coupled(:,:,c_data.P.coupled_select_ind,:)),...
            ones(T,1));
        B = c_data.B;
    catch
        sel_ind = selectParamCoupled(MM,P.theta_stds,top_dir,sim_name,alg_name,ps_name,nn);
        c_data = load(fullfile(top_dir,sim_name,alg_name,...
            [ps_name,'_',num2str(nn),'_',num2str(sel_ind),'.mat']));
        A0 = unshifted_basis_vector_ft_stack_zpad(c_data.P);
        fits = Ax_ft_1D_Time(A0,c_data.X_hat,ones(c_data.P.num_ims,1));
        load(fullfile(top_dir,sim_name,'IndepISM',...
             [ps_name,'_',num2str(nn),'.mat']),'B')
    end
    
    % axes(ha_waterfall(1))
    subplot(3,2,i)
    waterfall(B')
    zlim([0 1.5])
    ylabel('time')
    xlabel('\eta')
    zlabel('Intensity')
    i = i + 1;
    % axes(ha_waterfall(2))
    subplot(3,2,i)
    waterfall(fits(:,:)')
    zlim([0 1.5])
    ylabel('time')
    xlabel('\eta')
    zlabel('Intensity')
    i = i + 1;
end
figOut.Position = [80,80,400 500];

end


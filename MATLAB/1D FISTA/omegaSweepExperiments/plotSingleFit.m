% Look at particular fit
r = 3;
o = 4;
m = 40;
t = 31;
T = 546;
top_dir = 'E:\MMPAD_omega';
om_dir = {'omega2','omega3','omega4','omega5'};
r_dir = {'ring1','ring2','ring3','ring4'};


load(fullfile('C:\Users\dpqb1\Downloads',...
['indep_fit_',num2str(m),'_',num2str(t),'.mat']))

dset_name = r_dir{r};
om_name = om_dir{o};
load(fullfile(top_dir,om_name,dset_name,[P.prefix,'_',num2str(t),'.mat']));
bm = mirrorData(sum(polar_image,2));

A0 = unshifted_basis_vector_ft_stack_zpad(P);

% b = B(:,j);

figure(1)
fit = Ax_ft_1D(A0,x_hat);  
hold on
plot(bm*P.dataScale)
plot(fit)

az_signal = squeeze(sum(x_hat,1));
var_sum = sum(az_signal(:));
awmv = sum(sqrt(P.var_theta(:)).*az_signal(:))/var_sum

%% plot awmvs for original selected indices
awmv_dir = 'E:\MMPAD_omega\indep';
for r = 1:4
    figure(r)
    hold on
    for o = 1:4
        load(fullfile(awmv_dir,[r_dir{r},om_dir{o},'_mirror_indep_awmv']))
        awmv_select = zeros(T,1);
        for t = 1:T
            awmv_select(t) = awmv(t,select_indices(t));
        end
        plot(awmv_select)
    end
end

%% try different methods for parameter selection and plot awmvs
close all
top_dir = 'E:\MMPAD_omega';
om_dir = {'omega2','omega3','omega4','omega5'};
r_dir = {'ring1','ring2','ring3','ring4'};
awmv_dir = 'E:\MMPAD_omega\indep';
T = 546;
for r = 1:4
    figure(r)
    hold on
    % also plot awmvs from paper alongside
    p_data = load(fullfile('C:\Users\dpqb1\Desktop\AWMV_mirror_Figures',...
                [r_dir{r},'_zero_mirror_coupled_awmv']));
    c_data = load(fullfile('D:\MMPAD_data_nr1',...
        [r_dir{r},'_zero_coupled_CG_TVphi_Mirror5'],'coupled_fit_1'),'P');
    old_ind = c_data.P.params.lambda1_indices;
    
    plot(p_data.awmv_az_init,'Linewidth',2)
    hold on
    for o = 4
        load(fullfile(awmv_dir,[r_dir{r},om_dir{o},'_mirror_indep_awmv']))
        new_select_ind = zeros(T,1);
        for t = 1:T
            err_t = err_select(:,t);
            l1_t = l1_select(:,t);
            err_t = err_t;%/max(err_t);
            l1_t = l1_t;%/max(l1_t);
            sq_origin_dist = abs(l1_t) + abs(err_t);
            new_select_ind(t) = find(sq_origin_dist == min(sq_origin_dist));
        end
        
%         for t = 1:T
%             rel_err_t = rel_err_select(:,t);
%             while rel_err_t(new_select_ind(t)) > 0.2
%                 if new_select_ind(t) > 1
%                     new_select_ind(t) = new_select_ind(t) - 1;
%                 else
%                     new_select_ind(t) = find(rel_err_t == min(rel_err_t));
%                     break
%                 end
%             end
%         end
        
        
        awmv_select = zeros(T,1);
        awmv_old = zeros(T,1);
        exact_ind = zeros(T,1);
        for t = 1:T
            awmv_select(t) = awmv(t,new_select_ind(t));
            awmv_old(t) = awmv(t,old_ind(t));
%             exact_ind(t) = find(awmv(t,:) == p_data.awmv(t));
%             awmv_select(t) = awmv(t,exact_ind(t));
        end
        plot(awmv_select)
%         plot(awmv_old)
        
        
    end

    
    legend('paper','now','new-old ind')
end
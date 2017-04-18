function save_vars( dir,xhat_new,rel_fit_error_new,sparse_new,rel_fit_error_old,sparse_old )
%save_vars Summary of this function goes here
%   Detailed explanation goes here
save(dir,'xhat_new','rel_fit_error_new','sparse_new',...
                    'rel_fit_error_old','sparse_old')

end


%% Load fit data
outFile = 'spread_alphap_%2.2f_betap_%2.2f.mat';
result_path = fullfile('D:','CHESS_data','al7075_311_full_alphap_fit');
fileInfo = dir(fullfile(result_path,'fista_fit*.mat'));

for i = 596:596
    load(fullfile(fileInfo(i).folder,fileInfo(i).name))
    if(i == 1)
        var_signal = zeros(P.num_var_t,P.num_var_r,5,41,5);
        rel_error = zeros(5,41,5);
        sparsity = zeros(5,41,5);
    end
    disp(['Load: ' num2str(P.load_step) '   Image: ' num2str(P.img)])
    idx1 = floor(P.img/5)+1;
    idx2 = mod(P.img,5)+1;     
    sparsity(P.load_step+1,idx1,idx2) = sum(x_hat(:)>0);
    rel_error(P.load_step+1,idx1,idx2) = err(end);
    var_signal(:,:,P.load_step+1,idx1,idx2) = squeeze(sum(sum(x_hat,1),2));
end

save(sprintf(outFile,P.alphap,P.betap),'var_signal','rel_error','sparsity','P')
        
% Show simulated data

load("C:\Users\dpqb1\Documents\Dissertation\ms_sim2\output_j1_sig_1.00e-02_lam1_6.00e-02_lam2_0.00e+00.mat")
suffix = 'indep';
outputs_indep = outputs;

load("C:\Users\dpqb1\Documents\Dissertation\ms_sim2\output_j10_sig_1.00e-02_lam1_6.00e-02_lam2_5.00e-01.mat")
suffix = 'of';
outputs_of = outputs;
figDir = "C:\Users\dpqb1\Documents\Dissertation\ms_sim2";

Yhats_indep = fig2decomp(outputs_indep);
Yhats_of = fig2decomp(outputs_of);

Yhat1_indep = Yhats_indep{1};
Yhat2_indep = Yhats_indep{1};

Yhat1_of = Yhats_of{1};
Yhat2_of = Yhats_of{2};

[y,y_true,K,J,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center();
center = (M+1)/2;
Xtrue1 = Xtrue;
Xtrue2 = Xtrue;
Xtrue1(:,:,9:16,:) = 0;
Xtrue2(:,:,1:8,:) = 0;

scales = outputs.scales;

ADtrue = reSampleCustomArrayCenter(N,Dtrue,scales,center);
ADtrue = padarray(ADtrue,[0 M-1 0],0,'post');
ADftrue = fft2(ADtrue);

Ytrue1 = squeeze(ifft2(sum(bsxfun(@times,ADftrue,fft2(Xtrue1)),3),'symmetric'));
Ytrue1 = unpad(Ytrue1,M-1,'pre');

Ytrue2 = squeeze(ifft2(sum(bsxfun(@times,ADftrue,fft2(Xtrue2)),3),'symmetric'));
Ytrue2 = unpad(Ytrue2,M-1,'pre');

%%
for i = 1:30
    figure(12)
    subplot(2,1,1)
    plot(Ytrue1(:,i),'Linewidth',2)
    hold on
    plot(Ytrue2(:,i),'Linewidth',2)
    plot(Yhat1_indep(:,i),'Linewidth',2)
    plot(Yhat2_indep(:,i),'Linewidth',2)
    hold off


    subplot(2,1,2)
    plot(Ytrue1(:,i),'Linewidth',2)
    hold on
    plot(Ytrue2(:,i),'Linewidth',2)
    plot(Yhat1_of(:,i),'Linewidth',2)
    plot(Yhat2_of(:,i),'Linewidth',2)
    hold off
    pause()
end


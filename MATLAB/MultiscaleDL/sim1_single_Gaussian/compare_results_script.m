% Load true solution
sigmas = 0:0.01:0.1;
dataset = 'steps_matched';
[y,~,N,M,T,Xtrue,Dtrue] = gaus_example_switch_multiscale_dl(sigmas(1),dataset);

% Load results
load('C:\Users\dpqb1\Documents\Outputs2024_10_31__Dtrue1_Xzeros0\steps_matched_results_sig_1\output_j1_1_sig_0.00e+00_lam1_1.00e-04_lam2_0.00e+00.mat')
D = outputs.D;
X = outputs.X;
scales = outputs.scales;
N = outputs.N;
M = outputs.M;
y = squeeze(outputs.y);
K = outputs.K;
center = (M+1)/2;
Uarray = zeros(numel(scales),1);
for i = 1:numel(scales)
    Uarray(i) = size(scales{i},2);
end
Utotal = sum(Uarray);


Xf = fft2(X);
AD = reSampleCustomArrayCenter(N,D,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'));
Yhat = unpad(Yhat,M-1,'pre');

Jdf = 0.5*sum(vec(abs(Yhat-y).^2))/2;
Jl1 = sum(abs(X),'all');
Jobj = Jdf + Jl1;


Xf = fft2(Xtrue);
AD = reSampleCustomArrayCenter(N,Dtrue,scales,center);
AD = padarray(AD,[0 M-1 0],0,'post');
ADf = fft2(AD);
Yhat = squeeze(ifft2(sum(bsxfun(@times,ADf,Xf),3),'symmetric'));
Yhat = unpad(Yhat,M-1,'pre');

Jdf_true = 0.5*sum(vec(abs(Yhat-y).^2))/2;
Jl1_true  = sum(abs(Xtrue),'all');
Jobj_true  = Jdf_true  + Jl1_true ;



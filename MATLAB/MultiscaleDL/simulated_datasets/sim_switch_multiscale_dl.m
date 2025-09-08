function [yn,y,N,M,T,Xtrue,Dtrue] = sim_switch_multiscale_dl(sigma,dataset,plotFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    plotFlag = false;
end

switch dataset
    case 'unmatched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_unmatched_multiscale_dl(sigma);
    case 'matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_matched_multiscale_dl(sigma);
    case 'unmatched2'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_unmatched2_multiscale_dl(sigma);
    case 'steps'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps_multiscale_dl(sigma);
    case 'steps2'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps2_multiscale_dl(sigma);
    case 'steps_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps_matched_multiscale_dl(sigma);
    case 'steps_unmatched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps_unmatched_multiscale_dl(sigma);
    case 'peak_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = peak_example_steps_matched_multiscale_dl(sigma);
    case 'peak_matched2'
        [yn,y,N,M,T,Xtrue,Dtrue] = peak_example_steps_matched2_multiscale_dl(sigma);
    case 'sim2_gaussian_tooth_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = sim2_gaussian_tooth_matched(sigma);
    case 'sim2_gaussian_tooth_matched_mini'
        [yn,y,N,M,T,Xtrue,Dtrue] = sim2_gaussian_tooth_matched(sigma);
        T = 5;
        yn = yn(:,1:5);
        y = y(:,1:5);
        Xtrue = Xtrue(:,:,:,1:5);
    case 'sim2_tooth_backtooth_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = sim2_tooth_backtooth_matched(sigma);
    case 'sim2_gaussian_tooth_unmatched'
        [yn,y,N,M,T,Xtrue,Dtrue] = sim2_gaussian_tooth_unmatched(sigma);
    case 'sim2_gaussian_tooth_unmatched2'
        [yn,y,N,M,T,Xtrue,Dtrue] = sim2_gaussian_tooth_unmatched2(sigma);
    case 'sim2_gaussian_tooth_matched2'
        [yn,y,N,M,T,Xtrue,Dtrue] = sim2_gaussian_tooth_unmatched2(sigma);
    case 'dissertation'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center(sigma);
    case 'dissertation_adjust'
        [yn,y,N,M,T,Xtrue,Dtrue] = dissertation_adjust(sigma,plotFlag);
    case 'dissertation_adjust2'
        [yn,y,N,M,T,Xtrue,Dtrue] = dissertation_adjust2(sigma,plotFlag);
    case 'dissertation_adjust_nooverlap'
        [yn,y,N,M,T,Xtrue,Dtrue] = dissertation_adjust_nooverlap(sigma,plotFlag);
    case 'dissertation_long'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center_long(sigma,plotFlag);
    case 'dissertation_long_separate'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center_long_separate(sigma,plotFlag);
    case 'mmpad_ring1'
        ring_num = 1;
        try
            yn = loadMMPAD1D(ring_num,'interp','/cluster/home/dbanco02/mmpad');
        catch
            yn = loadMMPAD1D(ring_num,'interp','E:\Data\mmpad');
        end
        yn = yn(:,21:100);
        [N,T] = size(yn);
        y = yn;
        M = 181;
        Xtrue = 0;
        Dtrue = 0;
    case 'mmpad_ring3'
        ring_num = 3;
        try
            yn = loadMMPAD1D(ring_num,'copy-shift','/cluster/home/dbanco02/mmpad');
        catch
            yn = loadMMPAD1D(ring_num,'copy-shift','E:\Data\mmpad');
        end
        [N,T] = size(yn);
        y = yn;
        M = 51;
        Xtrue = 0;
        Dtrue = 0;

    case 'pseudo-voigt_unmatched'
        [yn,y,N,M,T,Xtrue,Dtrue] = voigt_example_unmatched_multiscale_dl(sigma);
    case 'pseudo-voigt_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = voigt_example_matched_multiscale_dl(sigma);
    case 'voigt_tooth_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = voigt_tooth_matched(sigma,plotFlag);
    case 'gaussian_tooth_matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaussian_tooth_matched(sigma,plotFlag);
    case 'gaussian_tooth_matched_long'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaussian_tooth_matched_long(sigma,plotFlag);
    case 'gaussian_tooth_matched_long2'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaussian_tooth_matched_long2(sigma,plotFlag);
    case 'voigt_tooth_matched_long3'
        [yn,y,N,M,T,Xtrue,Dtrue] = voigt_tooth_matched_long3(sigma,plotFlag);
    case 'voigt_example_6peaks'
        [yn,y,N,M,T,Xtrue,Dtrue] = voigt_example_6peaks(sigma,plotFlag);
    otherwise
        error('Invalid dataset name')
end

end
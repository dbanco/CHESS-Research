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
    case 'dissertation_long'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center_long(sigma,plotFlag);
    case 'dissertation_long_separate'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_linear_osc_signal_matched_small_zpad2_center_long_separate(sigma,plotFlag);
    case 'mmpad_ring1'
        ring_num = 1;
        yn = loadMMPAD1D(ring_num,1,'/cluster/home/dbanco02/mmpad');
        y = y(:,21:100);
        [N,T] = size(y);
        y = yn;
        M = N;
        Xtrue = 0;
        Dtrue = 0;
    otherwise
        error('Invalid dataset name')
end

end
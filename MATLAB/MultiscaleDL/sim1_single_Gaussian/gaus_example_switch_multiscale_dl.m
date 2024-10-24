function [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_switch_multiscale_dl(sigma,dataset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch dataset
    case 'unmatched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_unmatched_multiscale_dl(sigma);
    case 'matched'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_matched_multiscale_dl(sigma);
    case 'steps'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps_multiscale_dl(sigma);
    case 'steps2'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_steps2_multiscale_dl(sigma);
    case 'unmatched2'
        [yn,y,N,M,T,Xtrue,Dtrue] = gaus_example_unmatched2_multiscale_dl(sigma);
end

end
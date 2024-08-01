function D = initGaussianDictionary(dictFilterSizes)
K = size(dictFilterSizes,2);
P.Ny = max(dictFilterSizes(1,:));
P.Nx = max(dictFilterSizes(2,:));
P.stds = [logspace(log10(0.5),log10(P.Ny/8),K);...
          logspace(log10(0.5),log10(P.Nx/8),K)]';
P.means = [dictFilterSizes(1,:)/2;...
           dictFilterSizes(2,:)/2]'; 
D = dictionary2D(P);

end
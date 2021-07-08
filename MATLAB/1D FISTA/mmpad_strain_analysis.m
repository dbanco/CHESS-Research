top_dir = 'E:\MMPAD_data_nr1';
dset_name = 'ring%i_zero';
prefix = 'mmpad_img_%i.mat';

R = 4;
N = 546;
COM = zeros(R,N);
for ring = 1:R
    for img = 1:N
        ring_dir = sprintf(dset_name,ring);
        img_file = sprintf(prefix,img);
        load(fullfile(top_dir,ring_dir,img_file))
        
        [n,m] = size(polar_image);
        radial_img = sum(polar_image,2);
        radial_coord = (1:n)';
        total_mass = sum(radial_img);
        COM(ring,img) = sum(radial_img.*radial_coord./total_mass);
    end
end

plot(COM)
xlabel('time')
ylabel('center of mass (pixel)')
title('Radial Shifts in MMPAD CP-Ti')
legend('Ring 1','Ring 2','Ring 3','Ring 4','Location','Best')

figure(2)
plot(sum(COM))

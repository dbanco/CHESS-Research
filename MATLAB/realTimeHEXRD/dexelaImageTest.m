
% fname = '/nfs/chess/raw/2022-3/id1a3/miller-3528-a/ff-c103-90-s2-1/2/ff/ff1_000300.h5';
% r1 = 1195;
% r2 = 1225;
% center = [1944, 3560];
% [b, rad, az] = loadH5polar(1,1,fname,center,r1,r2);
% imagesc(b)

%ff1 vert flip
%ff2 hor filp
%10cm gap or 0.1m or 0.1e-6um

gap = round(10.5/0.0748);
center = [3026,1958];
fname = '/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/2/ff/ff1_000078.h5';
fname2 = '/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/2/ff/ff2_000078.h5';


h5disp(fname,'/imageseries/images')
img = squeeze(h5read( fname,'/imageseries/images', ...
                      [1 1 5],[3072 3888 40] ));
img2 = squeeze(h5read( fname2,'/imageseries/images', ...
                      [1 1 5],[3072 3888 40] ));

                 
% imagesc(full_dex) 
% 
% figure
% img_sum = squeeze(sum(img>500,[1,2]));
% [maxInt,maxInd] = max(img_sum);
% plot(img_sum)
% xlabel('\omega_i')
% title(sprintf('%0.2f at index %0.2f',maxInt,maxInd))

img_max = max(img,[],3);
img2_max = max(img2,[],3);
% full_dex = [img_max,zeros(3072,10),img2_max];
full_dex = [flip(img_max',1),zeros(3888,gap),flip(img2_max',2)]; 


figure(2)
imagesc(full_dex)
caxis([300 1000])
hold on

rad1 = 610*tan(pi*7.61/180)/0.0748 - 10;
rad2 = 610*tan(pi*7.61/180)/0.0748 - 55;
rad3 = 610*tan(pi*7.61/180)/0.0748 - 140;
for r = [rad1 rad2 rad3]
    r1 = r-20;
    r2 = r+20;
    

    theta = linspace(0,2*pi,300);

    x = r*cos(theta) + center(1);
    y = r*sin(theta) + center(2);
    x1 = r1*cos(theta) + center(1);
    y1 = r1*sin(theta) + center(2);
    x2 = r2*cos(theta) + center(1);
    y2 = r2*sin(theta) + center(2);
    
    plot(x,y,'y')
    plot(x1,y1,'g')
    plot(x2,y2,'g')
    
end
%%
% figure(1)
% imagesc(img)
% caxis([300 1000])
% hold on
% plot(x,y,'y')

% Show extracted rings in polar coordinates
figure(3)
i = 1;
for r = [rad1 rad2 rad3]
    r1 = r-20;
    r2 = r+20;
    [b, rad, az] = loadH5polarDex(14,14,fname,fname2,center,r1,r2);
    
    subplot(3,1,i)
    imagesc(b)
    caxis([300,2000])
    i = i + 1;
end



% hold on
% center = [1944,3560];
% r = 1215;
% theta = linspace(0,2*pi,300);
% 
% x = r*cos(theta) + center(1);
% y = r*sin(theta) + center(2);
% 
% plot(x,y)

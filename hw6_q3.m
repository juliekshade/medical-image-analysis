clear all
close all

im1 = imread('Slc1-06.gif');
im2 = imread('Slc1-07.gif');
im1 = double(im1);
im2 = double(im2);

n = 5;
d = 10;

SSD = nan(74,74);
u = nan(74,74);
v = nan(74,74);
pad = n+d;
im1p = padarray(im1, [pad pad], 0, 'both');
im2p = padarray(im2, [pad pad], 0, 'both');

for x = pad+1:pad+74
    for y = pad+1:pad+74
        SSD = zeros(2*d+1, 2*d+1);
        for dx = -d:d
            for dy = -d:d
                SSD_o = 0;
                for j = -n:n
                    for i = -n:n
                        SSD_o = SSD_o + ((double(im1p(x+i,y+j))-double(im2p(x+dx+i,y+dy+j)))^2);
                       
                    end 
                end
                SSD(dx+d+1,dy+d+1) = SSD_o;
             end    
        end
        [M,I] = min(SSD(:));
        [I_row, I_col] = ind2sub(size(SSD),I);
        u(x-pad,y-pad) = I_row - (d+1);
        v(x-pad,y-pad) = I_col - (d+1);
    end
end

[x,y] = meshgrid(1:74,1:74);
figure(3)
hold on;
quiver(x,y,u,v);
xlim([1 74]);
ylim([1 74]);
hax = gca; %get the axis handle
image(hax.XLim,hax.YLim,im1); %plot the image within the axis limits
quiver(x,y,u,v)
title('Velocity Field Solved with Block-Matching, n=5, d=10');

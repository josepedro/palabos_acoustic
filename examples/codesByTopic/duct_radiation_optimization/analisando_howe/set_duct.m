% colocar o duto nos caras

clear('all'); close all;clc;

load('min_ka/min_ka.mat');

acoustic_energy_duct = -acoustic_energy;
sizes = size(acoustic_energy);
times_all = sizes(3);
for time = 1:times_all
	acoustic_energy_duct(1:360, 24:25,time) = 0;
	acoustic_energy_duct(1:360, 58:59,time) = 0;
    acoustic_energy_duct(1:end, round(end/2):end,time) = fliplr(acoustic_energy_duct(1:end, 1:round(end/2),time));

    image_duct = acoustic_energy_duct(:, :,time);
    filter_image = fspecial('gaussian',[sizes(1) sizes(2)],0.8);
    image_duct = imfilter(image_duct,filter_image,'same');
    image_duct(1:360, 24:25) = 1;
    image_duct(1:360, 58:59) = 1;

	imagesc(flip(image_duct(:,:))), view(2), shading flat, caxis([-.0000001 .0000001]), axis equal
	grid off
	pause(.01) 

	% gif utilities
    set(gcf,'color','w'); % set figure background to white
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imiMd,cm] = rgb2ind(im,256);
    outfile = 'howe_min.gif';
 
    % On the first loop, create the file. In subsequent loops, append.
    if time==1
        imwrite(imiMd,cm,outfile,'gif','DelayTime',0.001,'loopcount',inf);
    else
        imwrite(imiMd,cm,outfile,'gif','DelayTime',0.001,'writemode','append');
    end
end

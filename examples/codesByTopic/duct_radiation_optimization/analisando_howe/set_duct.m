% colocar o duto nos caras

clear('all'); close all; clc;

load('max_ka_007/max_ka_007.mat');
load('max_ka_007/onda_acustica.mat');

acoustic_energy_duct = acoustic_energy;
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
    subplot(1,2,1);
    clims = [-.0000001 .0000001];
	imagesc(flip(image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside'), axis off
    %imagesc(flip(image_duct(320:end - 40,24:59))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside'), axis off
    txt3 = '2a/\pi \rightarrow';
    text(25,57,txt3,'HorizontalAlignment','right')
	%contour((image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside')
    grid off
    subplot(1,2,2);
    plot(b - mean(b), 'black'); hold on;
    plot(time, b(time) - mean(b), '*black')
	pause(.0000000001) 

	% gif utilities
    set(gcf,'color','w'); % set figure background to white
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imiMd,cm] = rgb2ind(im,256);
    outfile = 'max_ka_007.gif';
 
    % On the first loop, create the file. In subsequent loops, append.
    if time==1
        imwrite(imiMd,cm,outfile,'gif','DelayTime',0.001,'loopcount',inf);
    else
        imwrite(imiMd,cm,outfile,'gif','DelayTime',0.001,'writemode','append');
    end
end

acoustic_energy_duct_mean = sum(acoustic_energy_duct(329:end - 40,24:59, :),3)/sizes(3);

%acoustic_energy_duct_mean = sum(acoustic_energy_duct(280:end,:,:),3)/sizes(3);

acoustic_potent = sum(sum(acoustic_energy_duct_mean))

acoustic_energy_duct_mean = sum(acoustic_energy_duct(280:end,:,:),3)/sizes(3);
sizes_2 = size(acoustic_energy_duct_mean);
acoustic_energy_duct_mean(1:sizes_2(1) - 40, 24:25) = 1;
acoustic_energy_duct_mean(1:sizes_2(1) - 40, 58:59) = 1;

figure;imagesc(flip(acoustic_energy_duct_mean)), view(2), shading flat, caxis([-.0000001 .0000001]), axis equal, colormap(gray), colorbar('southoutside'), axis off
txt3 = '2a/\pi \rightarrow';
text(25,57,txt3,'HorizontalAlignment','right', 'FontSize', 25)

% figure;
% time = 190;
% %time = 402;
% %time = 232;
% %time = 46;
% %time = 108;
% %time = 72;
% %time = 37;
% image_duct = acoustic_energy_duct(:, :,time);
%     filter_image = fspecial('gaussian',[sizes(1) sizes(2)],0.8);
%     image_duct = imfilter(image_duct,filter_image,'same');
%     image_duct(1:360, 24:25) = 1;
%     image_duct(1:360, 58:59) = 1;
%     subplot(1,2,1);
%     clims = [-.0000001 .0000001];
% imagesc(flip(image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside'), axis off
%     txt3 = '2a/\pi \rightarrow';
%     text(25,57,txt3,'HorizontalAlignment','right','FontSize',20)
%     %contour((image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside')
%     grid off
%     subplot(1,2,2);
%     plot(b - mean(b), 'black'); hold on;
%     plot(time, b(time) - mean(b), 'oblack')
%     saveas(gcf,'max_007_1.png')


% figure;
% %time = 402;
% time = 300;
% %time = 232;
% %time = 46;
% %time = 108;
% %time = 72;
% %time = 37;
% image_duct = acoustic_energy_duct(:, :,time);
%     filter_image = fspecial('gaussian',[sizes(1) sizes(2)],0.8);
%     image_duct = imfilter(image_duct,filter_image,'same');
%     image_duct(1:360, 24:25) = 1;
%     image_duct(1:360, 58:59) = 1;
%     subplot(1,2,1);
%     clims = [-.0000001 .0000001];
% imagesc(flip(image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside'), axis off
%     txt3 = '2a/\pi \rightarrow';
%     text(25,57,txt3,'HorizontalAlignment','right','FontSize',20)
%     %contour((image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside')
%     grid off
%     subplot(1,2,2);
%     plot(b - mean(b), 'black'); hold on;
%     plot(time, b(time) - mean(b), 'oblack')
%     saveas(gcf,'max_007_2.png')

% figure;
% time = 409;
% %time = 402;
% %time = 232;
% %time = 46;
% %time = 108;
% %time = 72;
% %time = 37;
% image_duct = acoustic_energy_duct(:, :,time);
%     filter_image = fspecial('gaussian',[sizes(1) sizes(2)],0.8);
%     image_duct = imfilter(image_duct,filter_image,'same');
%     image_duct(1:360, 24:25) = 1;
%     image_duct(1:360, 58:59) = 1;
%     subplot(1,2,1);
%     clims = [-.0000001 .0000001];
% imagesc(flip(image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside'), axis off
%     txt3 = '2a/\pi \rightarrow';
%     text(25,57,txt3,'HorizontalAlignment','right','FontSize',20)
%     %contour((image_duct(280:end,:))), view(2), shading flat, caxis(clims), axis equal, colormap(gray), colorbar('southoutside')
%     grid off
%     subplot(1,2,2);
%     plot(b - mean(b), 'black'); hold on;
%     plot(time, b(time) - mean(b), 'oblack')
%     saveas(gcf,'max_007_3.png')

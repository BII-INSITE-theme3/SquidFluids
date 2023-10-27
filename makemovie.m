
%% Make a video 

v = VideoWriter("test.avi");
open(v)

% Loop through times
for i = 1:49

    filename = ['frames/flowAngle/' num2str(i) '.png'];
    A = imread(filename);
    imagesc(A)
    frame = getframe(gcf);
    writeVideo(v,frame);

    hold off

end

close(v)

%% Load the files

%filename = 'output_burger1.txt';
%filename = 'output_burger1.txt';
%%filename = ['output_firstOrder.txt'];
close all

for i = 0
filename = ['ds_convergence/test_' num2str(i) '.txt'];
A = readmatrix(filename);

plot(A(1,:))
hold on
end

%% Make a video 

filename = ['ds_convergence/test_' num2str(0) '.txt'];
A = readmatrix(filename);
% v = VideoWriter("secondOrderTest4.avi");
% open(v)

% Dimensionalisation
CFL = 0.5;
ds = 0.1;
dt = 0.005;

% Array sizes
[Nt,Ns] = size(A);

ss = linspace(0,ds*Ns,Ns); % Sizes
tt = linspace(0,dt*Nt,Nt); % Times

% Loop through times
for i = 1%:10:4000

    plot(ss,A(i,:),'k-','LineWidth',2)
    hold on;
    plot(ss,exp(-(ss-10-(i-1)*dt).^2),'r--','LineWidth',2)
    ylim([-1,10])

    % frame = getframe(gcf);
    % writeVideo(v,frame);

    hold off

end

% %% old
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v = VideoWriter("SPM.avi");
% open(v)

% range = 11:1:2000;
% s = length(A(1,:));
% c = jet(length(range));

% dt = 0.05;
% ds = dt*2;

% for i = 1:length(range) 

%     n = range(i);

%     for j = -10:10
%         plot(A(n+j,:),'color',c(i,:))
%         plot(sol())
%         ylim([0,1])
%         hold on;

%     end

%     frame = getframe(gcf);
%     writeVideo(v,frame);
%     hold off;

% end

% %%

% % Dimensionalisation
% CFL = 0.5;
% ds = 0.01;
% dt = ds*CFL;

% % Array sizes
% [Nt,Ns] = size(A);

% ss = linspace(0,ds*Ns,Ns); % Sizes
% tt = linspace(0,dt*Nt,Nt); % Times

% % Loop through times
% for i = 1:10:200

%     plot(ss,A(i,:),'k-','LineWidth',2)
%     hold on;
%     plot(ss,exp(-(ss-10-(i-1)*dt).^2),'r--','LineWidth',2)
%     hold off
%     pause(0.5)

% end
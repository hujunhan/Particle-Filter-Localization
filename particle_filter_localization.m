function [error] = particle_filter_localization(NP_out)
%PARTICLE_FILTER_LOCALIZATION Summary of this function goes here
%   Detailed explanation goes here

% -------------------------------------------------------------------------
% TASK for particle filter localization
% for robotic class in 2018 of ZJU

% Preparartion: 
% 1. you need to know how to code and debug in matlab
% 2. understand the theory of Monte Carlo

% Then complete the code by YOURSELF!
% -------------------------------------------------------------------------

%close all;
%clear all;

disp('Particle Filter program start!!')

%% initialization
time = 0;
endTime = 60; % second
global dt;
dt = 0.1; % second
 
nSteps = ceil((endTime - time)/dt);

localizer.time = [];
localizer.xEst = [];
localizer.xGnd = [];
localizer.xOdom = [];
localizer.z = [];
localizer.PEst=[];
localizer.u=[];
localizer.error=[];
% Estimated State [x y yaw]'
xEst=[0 0 0]';
% GroundTruth State
xGnd = xEst;
% Odometry-only = Dead Reckoning  
xOdom = xGnd;

% Covariance Matrix for predict
Q=diag([0.1 0.1 toRadian(3)]).^2;
% Covariance Matrix for observation
R=diag([1]).^2;% range:meter

% Simulation parameter
global Qsigma
Qsigma=diag([0.1 toRadian(5)]).^2;
global Rsigma
Rsigma=diag([0.1]).^2;

% landmark position
landMarks=[10 0; 10 10; 0 15; -5 20];


% longest observation confined
MAX_RANGE=20;
% Num of particles, initialized
NP=NP_out;
% Used in Resampling Step, a threshold
NTh=NP/2.0;

% particles produced
px=repmat(xEst,1,NP);
% weights of particles produced
pw=zeros(1,NP)+1/NP;


%% Main Loop 

    
for i=1 : nSteps
    %disp(size(px))
    time = time + dt;
    u=doControl(time);
    
    % do observation
    [z,xGnd,xOdom,u]=doObservation(xGnd, xOdom, u, landMarks, MAX_RANGE);
    error=0;
    for ip=1:NP
        
        % process every particle
        x=px(:,ip);
        %w=1;
        w=pw(:,ip);
        dx=x(1)-xGnd(1);%Îó²î¼ÆËã
        dy=x(2)-xGnd(2);
        error=error+sqrt(dx^2+dy^2);
        % do motion model and random sampling
        x=doMotion(x, u)+sqrt(Q)*randn(3,1);
         % calculate inportance weight
        for iz=1:length(z(:,1))
            pz=norm(x(1:2)'-z(iz,2:3));
            dz=pz-z(iz,1);
            w=w*Gaussian(dz,0,sqrt(R));
        end
        px(:,ip)=x;
        pw(ip)=w;
        
    end
    error=error/NP;
    pw=Normalization(pw,NP);
    xEst=px*pw';
    [px,pw]=ResamplingStep(px,pw,NTh,NP);

    
    
    % Simulation Result
    localizer.time=[localizer.time; time];
    localizer.xGnd=[localizer.xGnd; xGnd'];
    localizer.xOdom=[localizer.xOdom; xOdom'];
    localizer.xEst=[localizer.xEst;xEst'];
    localizer.u=[localizer.u; u'];
     localizer.error=[localizer.error,error];
    %Animation (remove some flames)
    subplot(2,1,1)
    if rem(i,10)==0 
        hold off;
        arrow=0.5;
        for ip=1:NP
            quiver(px(1,ip),px(2,ip),arrow*cos(px(3,ip)),arrow*sin(px(3,ip)),'ok');hold on;
        end
        plot(localizer.xGnd(:,1),localizer.xGnd(:,2),'.b');hold on;
        plot(landMarks(:,1),landMarks(:,2),'pk','MarkerSize',10);hold on;
        if~isempty(z)
            for iz=1:length(z(:,1))
                ray=[xGnd(1:2)';z(iz,2:3)];
                plot(ray(:,1),ray(:,2),'-r');hold on;
            end
        end
        plot(localizer.xOdom(:,1),localizer.xOdom(:,2),'.k');hold on;
        plot(localizer.xEst(:,1),localizer.xEst(:,2),'.r');hold on;
        axis equal;
        grid on;
        drawnow;
    end
    
end
error=localizer.error;
N=length(error);
error=sum(error)/N;
% 
%draw the final results of localizer, compared to odometry & ground truth
drawResults(localizer);
subplot(2,1,2);
plot(localizer.error)
title('Average Error', 'fontsize', 12, 'fontname', 'times');
xlabel('Time(dt)', 'fontsize', 12, 'fontname', 'times');
ylabel('Average Error (m)', 'fontsize', 12, 'fontname', 'times');
end









%% Other functions

% degree to radian
function radian = toRadian(degree)
    radian = degree/180*pi;
end

function []=drawResults(localizer)
%Plot Result
 
    figure(1);
    hold off;
    x=[ localizer.xGnd(:,1:2) localizer.xEst(:,1:2)];
    set(gca, 'fontsize', 12, 'fontname', 'times');
    plot(x(:,1), x(:,2),'-.b','linewidth', 4); hold on;
    plot(x(:,3), x(:,4),'r','linewidth', 4); hold on;
    plot(localizer.xOdom(:,1), localizer.xOdom(:,2),'--k','linewidth', 4); hold on;

    title('Localization Result', 'fontsize', 12, 'fontname', 'times');
    xlabel('X (m)', 'fontsize', 12, 'fontname', 'times');
    ylabel('Y (m)', 'fontsize', 12, 'fontname', 'times');
    legend('Ground Truth','Particle Filter','Odometry Only');
    grid on;
    axis equal;

end

function [ u ] = doControl( time )
%DOCONTROL Summary of this function goes here
%   Detailed explanation goes here

    %Calc Input Parameter
    T=10; % [sec]

    % [V yawrate]
    V=1.0; % [m/s]
    yawrate = 5; % [deg/s]

    u =[ V*(1-exp(-time/T)) toRadian(yawrate)*(1-exp(-time/T))]';


end


%%  you need to complete

% do Observation model 
function [z, xGnd, xOdom, u] = doObservation(xGnd, xOdom, u, landMarks, MAX_RANGE)
    global Qsigma;
    global Rsigma;
    
    % Gnd Truth and Odometry
    xGnd=doMotion(xGnd, u);% Ground Truth ÀíÏë×´Ì¬
    u=u+sqrt(Qsigma)*randn(2,1);% add noise randomly
    xOdom=doMotion(xOdom, u); % odometry only
    
    %Simulate Observation
    z=[];
    for iz=1:length(landMarks(:,1))
        dx = xGnd(1)-landMarks(iz,1);
        dy = xGnd(2)-landMarks(iz,2);
        d=sqrt(dx^2+dy^2);
        if d<MAX_RANGE 
            z=[z;[d+sqrt(Rsigma)*randn(1,1) landMarks(iz,:)]];   % add observation noise randomly
        end
    end
end


% do Motion Model
function [ x ] = doMotion( x, u)
    global dt;
    
    Delta = [ [dt*cos(x(3)+u(2)),0];
              [dt*sin(x(3)+u(2)),0];
              [0,dt]];

    x = x+Delta*u;
end

% Gauss function
function g = Gaussian(x,u,sigma)
    g=exp(-((x-u)^2)/((sigma^2)*2.0))/sqrt(2.0*pi*(sigma^2));
end

% Normalization 
function pw=Normalization(pw,NP)
    pw=pw/sum(pw);

end

% Resampling
function [px,pw]=ResamplingStep(px,pw,NTh,NP)
    Neff=1.0/(pw*pw');
    %Neff=0;
    if Neff<NTh
        ww=pw(1);
        for iw=2:NP
            ww=[ww,ww(end)+pw(iw)];
        end
        pw1=[];
        pp=[];
        for i=1:NP
            r=rand();
            for j=1:NP
                if ww(j)>r
                    pp=[pp,px(:,j)]; 
                    pw1=[pw1,pw(:,j)];
                    break
                end
            end
        end
        px=pp;
        pw=pw1;
    end
end
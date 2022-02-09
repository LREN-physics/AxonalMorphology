function EstimateMicrostructure_GroupDelay(subject_list,path,path_group_results,delay,ci)
% Estimates morphological features of white matter tracts. In particular,
% the axon radius distribution, P(r), and the g-ratio dependence on the
% radius, g(r).
% P(r) is a gamma distribution defined by two parameters: theta and mode
% (both in micrometres). Theta is estimated here from the data. M is fixed
% to 0.4.
% G(r) is a power-law defined by two parameters: alpha (unitless) and beta
% (micrometres^-alpha). Alpha is fixed to 0.14 and beta is estimated from
% the data here.
%
%   INPUTS:
%       subject_list (cell array)       - list of subject we want to run
%       path (string)                   - path to the subjects data folder
%       path_group_results (string)     - path to group results folder
%       delay (number)                  - IHTT in miliseconds
%       ci (number)                     - confidence interval in IHTT
%         estimation obtained in milisecond
%
%   OUPUTS:
%       none
%
% Author: Rita Oliveira
% Email: rita.oliveira.uni@gmail.com
% Laboratory for Research in Neuroimaging (LREN)
% Department of Clinical Neuroscience, Lausanne  University Hospital and University of Lausanne
% Mont-Paisible 16, CH-1011 Lausanne, Switzerland
%
% Last updated: 11/11/2021
%
%------------- BEGIN CODE --------------

close all

% Define colors and open figures and text files
Color_list=[230, 25, 75;60, 180, 75;255, 225, 25;0, 130, 200;245, 130, 48; 70, 240, 240;240, 50, 230;250, 190, 212;0, 128, 128;220, 190, 255;170, 110, 40;255, 250, 200;128, 0, 0;170, 255, 195;0, 0, 128;128, 128, 128;0, 0, 0]./255;
f3=figure(3);
set(f3,'units','normalized','outerposition',[0 0 1 1]);
ax3= axes('Parent',f3);
set(f3','visible','off');
    
% Initiate variables
AxonalDistribution_FixM_Group = [];

% LOOP FOR SUBJECTS
for sub=1:length(subject_list)
    
    % Define initial stuff
    sub_name = subject_list{sub};
    sub_dir  = fullfile(path,sub_name,'MRI');
    fprintf('Working on subject %s \n',sub_name);
    AxonalDistribution_FixM = [];
    
    % Start new figure
    if sub~=1
        close(f1)
        close(f2)
    end
    f1=figure(1);
    set(f1,'units','normalized','outerposition',[0 0 1 1]);   
    ax1= axes('Parent',f1);
    set(f1','visible','off');
    
    f2=figure(2);
    set(f2,'units','normalized','outerposition',[0 0 1 1]);   
    set(f2','visible','off');
    ax2= axes('Parent',f2);
    
    % Save information
    if sub==1
        AxonalDistribution_FixM_Group.SubNames            = {sub_name};
    else
        AxonalDistribution_FixM_Group.SubNames{end+1}     = sub_name;
    end
    
    % Create folders
    path_model          = fullfile(sub_dir, 'MicrostructureEstimation');
    if~exist(path_model,'dir')
        mkdir(path_model)
    end    
   
    % Load files for G-Ratio and tract length
    load(fullfile(sub_dir,'G_ratio_samples.mat'),'G');
    load(fullfile(sub_dir,'Tract_length.mat'),'len');
    
    % Velocity
    v = len/delay;
    save(fullfile(sub_dir,'Velocity.mat'),'v')
    
    % Initiate
    alpha   = 0.14;         % alpha from g(r)
    M       = 0.4;          % M from P(r), i.e, mode/peak of the distribution
    p       = 5.5;          % proportionality factor that relates the fiber diameter (D) with the conduction velocity in: V=p.D
    
    % Options for solver
    x0      = [0.7;0.1]; % First guess on M and theta
    options = optimoptions('lsqnonlin', 'Display','off','FunctionTolerance',10^-14,'OptimalityTolerance',10^-14,'algorithm','trust-region-reflective');
    
    %% >> Get model results using the limits of the confidence interval

    % Initiate
    clear Theta_interval Beta_interval p_r g_r
    X     = (0:0.05:5);                         % range of possible radius (micrometers)
    v_vec = len/(delay+ci):0.05:len/(delay-ci);  % range of velocities according to the std of the delay obtained
    if isempty(v_vec)
        disp('PROBLEM: no velocity interval. Check values for confidence interval and IHTT')
    end
    
    for j=1:length(v_vec)
        
        % Solve
        res     = lsqnonlin(@(x) MODEL_IN_alpha_mode_OUT_beta_theta(x,alpha,M,G,v_vec(j),p),x0,[0 0],[1 1],options); % other option
        beta        = res(1);
        theta       = res(2);
        
        % Get results
        shape       = M/theta+1;
        p_r         = gampdf(X,shape,theta);
        g_r         = beta.*X.^alpha;
        
        % Plot p(r)
        plot(ax1,X,p_r,'Color',[207 207 207]/250,'Linewidth',2); 
        hold(ax1,'on')
        
        % Plot g(r)
        plot(ax2,X,g_r,'Color',[207 207 207]/250,'LineWidth',2);
        hold(ax2,'on')
     
        % Save interval      
        if j==1 
            Velocity_interval(1) = v_vec(1);
            Theta_interval(1)    = theta;
            Beta_interval(2)     = beta;
        elseif j==length(v_vec)
            Velocity_interval(2) = v_vec(length(v_vec));
            Theta_interval(2)    = theta;
            Beta_interval(1)     = beta;
        end
        
    end

    %% >> Get model results using the main IHTT value
    
    % Solve
    res     = lsqnonlin(@(x) MODEL_IN_alpha_mode_OUT_beta_theta(x,alpha,M,G,v,p),x0,[0 0],[1 1],options); % other option
    beta        = res(1);
    theta       = res(2);
    
    % Get results
    shape       = M/theta+1;
    p_r         = gampdf(X,shape,theta);
    g_r         = beta.*X.^alpha;
    
    % Print results
    fprintf('>> Visual Transcallosal tract: \n Beta is: %0.02f um^-alpha \n Theta is: %0.02f um \n Velocity is: %0.02f m/s \n',beta,theta,v);
    
    % Plot main p(r)
    plot(ax1,X,p_r,'Color','k','Linewidth',2);
    ylim(ax1,[0 1.4])
    xlim(ax1,[0 5])
    xlabel(ax1,'Axon radius ({\mu}m)')
    ylabel(ax1,'P(r)')
    set(ax1,'box','off','color','none')
    set(ax1,'FontName','Helvetica','Fontsize',12)
    saveas(f1,fullfile(path_model,'P_r.png'))
     
    % Plot main g(r)
    plot(ax2,X,g_r,'Color','k','LineWidth',2);
    xlabel(ax2,'Axon radius ({\mu}m)')
    ylabel(ax2,'g(r)')
    set(ax2,'box','off','color','none')
    set(ax2,'FontName','Helvetica','Fontsize',12)
    saveas(f2,fullfile(path_model,'G_r.png'))

    % Plot g(r) with all the subjects together
    %figure(3)
    vec(sub)=plot(ax3,X,g_r,'Color',Color_list(sub,:),'LineWidth',2);
    hold(ax3,'on')
    
    % Save results
    AxonalDistribution_FixM.Beta              = beta;
    AxonalDistribution_FixM.Theta             = theta;
    AxonalDistribution_FixM.Velocity          = v;
    AxonalDistribution_FixM.Beta_interval     = Beta_interval;
    AxonalDistribution_FixM.Theta_interval    = Theta_interval;
    AxonalDistribution_FixM.Velocity_interval = Velocity_interval;
    AxonalDistribution_FixM.TractName         = 'V1V2';
    AxonalDistribution_FixM.M                 = M;
    AxonalDistribution_FixM.Alpha             = alpha;
    AxonalDistribution_FixM.P_factor          = p;
    save(fullfile(path_model,'AxonalDistribution_FixM.mat'),'AxonalDistribution_FixM')

end % end subjects

%% >> Finish

% Save image of g(r) for all subjects
figure(3)
xlabel('Axon radius ({\mu}m)','FontSize',12)
ylabel('g(r)','FontSize',12)
title('g(r) for all subjects','FontSize',14)
ylim([0.4 0.9])
legend(vec,subject_list,'Interpreter','None','Location','SouthEast')
saveas(gcf,fullfile(path_group_results,'G_r_group.png'))




end


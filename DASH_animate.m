% This script animates a given DASH leg module. All parameters are provided
% in as an input state vector x. v includes the animation video details and
% a boolean variable requiring a animation video output if true.

% % Each x element is shown below:
% % x(1)  -- l_1                               % link lengths of lift
% % x(2)  -- l_2
% % x(3)  -- l_3
% % x(4)  -- l_4                               % link lengths of swing
% % x(5)  -- l_5
% % x(6)  -- l_6
% % x(7)  -- theta_1                           % cum joint angles of lift
% % x(8)  -- theta_1 + theta_2
% % x(9)  -- theta_1 + theta_2 + theta_3
% % x(10) -- theta_4                           % cum joint angles of swing
% % x(11) -- theta_4 + theta_5
% % x(12) -- theta_4 + theta_5 + theta_6
% % x(13) -- L                                 % Length of the leg
% % x(14) -- l                                 % Length of coupled chain

% Function call
function out = DASH_animate(x, v, p)

% Unpack the kinematic state variables
l1 = x(1); 
l2 = x(2); 
l3 = x(3); 
l4 = x(4); 
l5 = x(5); 
l6 = x(6);
theta1 = x(7); 
theta2 = x(8) - theta1; 
theta3 = x(9) - theta2 - theta1;
theta4 = x(10); 
theta5 = x(11) - theta4; 
theta6 = x(12) - theta4 - theta5;
legL = x(13); kinL = x(14);

% ime vector for the simulation
t = linspace(0,1); out.t = t;

% % % LEGACY CODE, NOT SUPPORTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Conserved motor inputs to the system, scaled based on the length of the
% % % coupled kinematic chain, x(14) or kinL. The motor input is based on the
% % % initial condition given the optimization that has a kinL of 1. If the
% % % isolate input flag isoF is set to 1 (or 2) we only set the lift input (or
% % % swing input).
% % if isnan(iF) % default circular motor input case
% %     b_hat = deg2rad(15)/kinL;
% %     a_hat = b_hat;
% % elseif iF == 1 % just lift input case
% %     b_hat = 0;
% %     a_hat = deg2rad(15)/kinL;
% % else % swing input case
% %     b_hat = deg2rad(15)/kinL;
% %     a_hat = 0;
% % end
% The motor inputs to the system
b_hat = deg2rad(15)/kinL;
a_hat = b_hat;


% Load our analytical kinematics results -- manually solved kinematics
load('DASH_analytical.mat');

% Create symbols to help interpret the analytical functions
syms l_1 l_2 l_3 l_4 l_5 l_6 l L real positive
syms theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 real
syms tau real
syms alpha_hat beta_hat real positive

% Get each frame description to plot the kinematics ~~~~~~~~~~~~~~~~~~~~~~~
% Setup the error objective choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p struct incorporates this; p.io == 0 --> simpler formulation
%                             p.io == 1 --> complex formulation
%                                   p.c == 0 --> doubly coupled
%                                   p.c == 1 --> S lever, L coupled
%                                   p.c == 2 --> S coupled, L lever
h_e__i_f = matlabFunction(h_e__i, 'Vars', [tau, alpha_hat, beta_hat, l]);
h_i__1_t = double(subs(h_i__1, theta_1, theta1));
h_i__2_t = double(subs(h_i__1*h_1__2, [theta_1, l_1, theta_2], [theta1, l1, theta2]));
h_i__3_t = double(subs(h_i__1*h_1__2*h_2__3, [theta_1, l_1, theta_2, l_2, theta_3], [theta1, l1, theta2, l2, theta3]));
h_i__0bl_t = double(subs(h_i__1*h_1__2*h_2__3*h_3__0b_l, [theta_1, l_1, theta_2, l_2, theta_3, l_3], [theta1, l1, theta2, l2, theta3, l3]));
h_i__4_t = double(subs(h_i__4, theta_4, theta4));
h_i__5_t = double(subs(h_i__4*h_4__5, [theta_4, l_4, theta_5], [theta4, l4, theta5]));
h_i__6_t = double(subs(h_i__4*h_4__5*h_5__6, [theta_4, l_4, theta_5, l_5, theta_6], [theta4, l4, theta5, l5, theta6]));
h_i__0bs_t = double(subs(h_i__4*h_4__5*h_5__6*h_6__0b_s, [theta_4, l_4, theta_5, l_5, theta_6, l_6], [theta4, l4, theta5, l5, theta6, l6]));
h_i__23_t = double(subs(h_i__23, [theta_1, l_1, theta_2, l_2], [theta1, l1, theta2, l2]));
h_i__56_t = double(subs(h_i__56, [theta_4, l_4, theta_5, l_5], [theta4, l4, theta5, l5]));
switch p.c % Get the right transforms based on the coupling condition.
    case 0
        h_i__oprime_t = double(subs(h_i__oprime, [theta_1, l_1, theta_2, l_2, theta_4, l_4, theta_5, l_5], [theta1, l1, theta2, l2, theta4, l4, theta5, l5]));
        h_i__o_t = double(subs(h_i__o, [theta_1, l_1, theta_2, l_2, theta_4, l_4, theta_5, l_5, L], [theta1, l1, theta2, l2, theta4, l4, theta5, l5, legL]));
    case 1
        h_i__oprime_t = double(subs(h_i__oprimeL, [theta_1, l_1, theta_2, l_2, theta_4, l_4, theta_5, l_5], [theta1, l1, theta2, l2, theta4, l4, theta5, l5]));
        h_i__o_t = double(subs(h_i__oL, [theta_1, l_1, theta_2, l_2, theta_4, l_4, theta_5, l_5, L], [theta1, l1, theta2, l2, theta4, l4, theta5, l5, legL]));
    case 2
        h_i__oprime_t = double(subs(h_i__oprimeS, [theta_1, l_1, theta_2, l_2, theta_4, l_4, theta_5, l_5], [theta1, l1, theta2, l2, theta4, l4, theta5, l5]));
        h_i__o_t = double(subs(h_i__oS, [theta_1, l_1, theta_2, l_2, theta_4, l_4, theta_5, l_5, L], [theta1, l1, theta2, l2, theta4, l4, theta5, l5, legL]));
end

% LEGACY CODE, not supported anymore. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Compute the coupled length of the kinematic chain -- the average between
% % % % the swing and lift chains' net length (they should be close to equal due
% % % % to the equality constraints we have imposed on the optimizer, but it is
% % % % better to average)
% % % sRad = (norm(h_i__0bl_t(1:3,4))+norm(h_i__0bs_t(1:3,4)))*0.5;

% All frames are plotted using quiver. 
% 
% % Create empty containers to store the frame data as a function of time.
% The transforms below just shown the final frame since the starting from
% is the group origin/ identity e (dropped for convenience)
h_i_t = zeros(4,4,length(t));
h_1_t = h_i_t;
h_2_t = h_i_t;
h_3_t = h_i_t;
h_4_t = h_i_t;
h_5_t = h_i_t;
h_6_t = h_i_t;
h_0b_t = h_i_t;
h_23_t = h_i_t;
h_56_t = h_i_t;
h_oprime_t = h_i_t;
h_o_t = h_i_t;

% Compute the input frame location before running the video loop
for i = 1:length(t)
    % Obtain the location of the input frame; use that to compute
    % everything else
    h_i_t(:,:,i) = h_e__i_f(t(i), a_hat, b_hat, kinL);
    h_1_t(:,:,i) = h_i_t(:,:,i)*h_i__1_t; % lift frames
    h_2_t(:,:,i) = h_i_t(:,:,i)*h_i__2_t;
    h_3_t(:,:,i) = h_i_t(:,:,i)*h_i__3_t;
    h_4_t(:,:,i) = h_i_t(:,:,i)*h_i__4_t; % swing frames
    h_5_t(:,:,i) = h_i_t(:,:,i)*h_i__5_t;
    h_6_t(:,:,i) = h_i_t(:,:,i)*h_i__6_t;
    h_0b_t(:,:,i) = h_i_t(:,:,i)*(h_i__0bs_t + h_i__0bl_t)*0.5; % averaged from the two
    h_23_t(:,:,i) = h_i_t(:,:,i)*h_i__23_t; % mid-lift chain
    h_56_t(:,:,i) = h_i_t(:,:,i)*h_i__56_t; % mid-swing chain
    h_oprime_t(:,:,i) = h_i_t(:,:,i)*h_i__oprime_t; % leg base
    h_o_t(:,:,i) = h_i_t(:,:,i)*h_i__o_t; % leg-tip
end

% Return the input and output trajectory
out.h_i_t = h_i_t;
out.h_o_t = h_o_t;

% If a video of the animation is needed, then proceed.
if v.vidF
    % Colors for each chain
    circ1 = (1/255)*[215, 25, 28]; % input
    circ3 = (1/255)*[44, 123, 182]; % output
    circ2 = (1/255)*[96, 0, 220]; % lift
    circ4 = (1/255)*[35,139,69]; % swing
    
    % Plotting parameters (DEFINE BELOW)
    qW_f = 2.00; % coordinate frames line width
    qW_c = 2.50; % kinematic chain line width
    qW_L = 2.50; % Leg line width
    qW_t = 3.00; % input and output trajectory line width
    arrwL = 0.20; % arrow length scaling
    revS = 50; % try filled or unfilled for the 2D 
    fontS = 17; % font size for the assisting text
    
    % Let's rotate the camera a complete 360degs:
    azi_cosine = -37.5*cos(2*pi*t); % cosinusoidally oscillate the side view of camera
    % azi_sine = -37.5*sin(2*pi*t); % sinusoidally oscillate
    % azi_linear = 360*t; % a complete rotation
    
    % Create the video object
    video = VideoWriter(strcat(v.name,'.mp4'),'MPEG-4');
    % Set the framerate
    video.FrameRate = round(length(t)/v.vidT); % roughly takes about 10seconds to complete the video
    % Set the quality
    video.Quality = v.Quality;
    % Open the video
    open(video);
    
    % Start plotting
    figure('units','pixels','position',[0 0 1920 1080],'Color','w')
    
    % % Set the background color as slightly gray for maximum contrast between
    % % various elements in the figure
    % set(gcf, 'color', (1/255)*[210,210,210]);
    
    % Iterate over each frame
    for i = 1:length(t)
    
        % Not needed right now -- doesn't add anything to the plot
    % %     % Plot the origin
    % %     quiver3(h_e_t(1,4)*ones(1,3), h_e_t(3,4)*ones(1,3), h_e_t(2,4)*ones(1,3),...
    % %         arrwL*h_e_t(1,1:3), arrwL*h_e_t(3,1:3), arrwL*h_e_t(2,1:3), 0,...
    % %         'LineWidth', qW_f, 'Color', 'k'); % x-axis at the origin
    
        % Axis settings
        hold on; axis equal; axis([-0.5, 4.75, -2.0, 2.0, -1.5, 4.75]); %[-0.5, 4.75, -1.5, 1.5, -0.5, 4.75] % old limit
        axis off;
    %     grid on;
    
        % Plot the input %%%%%%%%% INPUT
        quiver3(h_i_t(1,4,i)*ones(1,3), h_i_t(3,4,i)*ones(1,3), h_i_t(2,4,i)*ones(1,3),...
            arrwL*h_i_t(1,1:3,i), arrwL*h_i_t(3,1:3,i), arrwL*h_i_t(2,1:3,i), 0,...
            'LineWidth', qW_f, 'Color', circ1); % coordinate axis at input
    
        % Plot the lift and swing kinematic chains
        plot3([h_i_t(1,4,i) h_1_t(1,4,i) h_2_t(1,4,i) h_3_t(1,4,i) h_0b_t(1,4,i)],...
            [h_i_t(3,4,i) h_1_t(3,4,i) h_2_t(3,4,i) h_3_t(3,4,i) h_0b_t(3,4,i)],...
            [h_i_t(2,4,i) h_1_t(2,4,i) h_2_t(2,4,i) h_3_t(2,4,i) h_0b_t(2,4,i)],...
            'LineWidth', qW_c, 'Color', circ2); % lift chain % , 'Marker', '^'
        plot3([h_23_t(1,4,i) h_oprime_t(1,4,i)],...
            [h_23_t(3,4,i) h_oprime_t(3,4,i)],...
            [h_23_t(2,4,i) h_oprime_t(2,4,i)],...
            'LineWidth', qW_c, 'Color', circ2); % lift chain to leg connection
        plot3([h_i_t(1,4,i) h_4_t(1,4,i) h_5_t(1,4,i) h_6_t(1,4,i) h_0b_t(1,4,i)],...
            [h_i_t(3,4,i) h_4_t(3,4,i) h_5_t(3,4,i) h_6_t(3,4,i) h_0b_t(3,4,i)],...
            [h_i_t(2,4,i) h_4_t(2,4,i) h_5_t(2,4,i) h_6_t(2,4,i) h_0b_t(2,4,i)],...
            'LineWidth', qW_c, 'Color', circ4); % swing chain % , 'Marker', '>'
        plot3([h_56_t(1,4,i) h_oprime_t(1,4,i)],...
            [h_56_t(3,4,i) h_oprime_t(3,4,i)],...
            [h_56_t(2,4,i) h_oprime_t(2,4,i)],...
            'LineWidth', qW_c, 'Color', circ4); % swing chain to leg connection
    
        % Plot the mechanical ground frame
        quiver3(h_0b_t(1,4,i)*ones(1,3), h_0b_t(3,4,i)*ones(1,3), h_0b_t(2,4,i)*ones(1,3),...
            1.5*arrwL*h_0b_t(1,1:3,i), 1.5*arrwL*h_0b_t(3,1:3,i),...
            1.5*arrwL*h_0b_t(2,1:3,i), 0,...
            'LineWidth', qW_f, 'Color', 'k'); % coordinate axis at gnd
    
        % Scatter a circle at mechanical ground (2 DOF revolute joint)
        scatter3(h_0b_t(1,4,i), h_0b_t(3,4,i), h_0b_t(2,4,i), revS,...
            'k', 'LineWidth', qW_f); %, 'filled'
        % unfilled circle is better to showcase a revolute joint
    
        % Plot the leg-base and leg tips using quiver %%%%%%%%% OUTPUT
        quiver3(h_oprime_t(1,4,i), h_oprime_t(3,4,i), h_oprime_t(2,4,i),...
            h_o_t(1,4,i)-h_oprime_t(1,4,i), h_o_t(3,4,i)-h_oprime_t(3,4,i),...
            h_o_t(2,4,i)-h_oprime_t(2,4,i), 0,...
            'LineWidth', qW_L, 'Color', circ3); % coordinate axis at gnd
    
        % Plot the trajectory here:
        plot3(reshape(h_i_t(1,4,:),1,length(t)),...
            reshape(h_i_t(3,4,:),1,length(t)),...
            reshape(h_i_t(2,4,:),1,length(t)),...
            'Color', circ1, 'LineWidth', qW_t); % input traj
        plot3(reshape(h_o_t(1,4,:),1,length(t)),...
            reshape(h_o_t(3,4,:),1,length(t)),...
            reshape(h_o_t(2,4,:),1,length(t)),...
            'Color', circ3, 'LineWidth', qW_t); % output traj
        
        % Labeling and titling
    %     xlabel('x'); ylabel('y'); zlabel('z');
    %     title('DASH 3D Kinematics');
    
        % Add text to denote different things
        text(-0.75,-0.25,-0.5,'Motor Input','Interpreter','latex',...
            'FontSize',fontS+5,'Color',circ1,...
            'FontWeight','bold')
    %     text(1.5,0.0,0.5,'Lift Chain','Interpreter','latex',...
    %         'FontSize',fontS-2,'Color',circ2)
        text(0.5,1.5,1,'Leg','Interpreter','latex',...
            'FontSize',fontS-2,'Color',circ3)
        text(4,1.5,2.25,'Output','Interpreter','latex',...
            'FontSize',fontS+5,'Color',circ3,...
            'FontWeight','bold')
    %     text(-0.65,1.65,-0.25,'Swing Chain','Interpreter','latex',...
    %         'FontSize',fontS-2,'Color',circ4)
    %     text(1.6,0.15,-0.5,'Ground','Interpreter','latex',...
    %         'FontSize',fontS-2,'Color','k')
        
        % Change the view: %azi_linear(i)
        view(azi_cosine(i),30); % sinusoidal camera rotation
        % azi_sine % azi_linear
    
    %     view(azi_linear(i),30); % linear full 360degs camera rotation
    %     view(-45,30); % default values -- good starting point
    %     view(-37.5,135); % beyond top view
    %     view(-90,30); % top-ish view
    %     view(0,90); % top view - SWING
    %     view(0,0); % front view - LIFT
    
        % Draw the current plot since we are in a loop:
        drawnow;
        % Capture the frame:
        writeVideo(video,getframe(gcf));
        % Clear the frame if it is not the final one:
        if i ~= length(t)
            clf;
        end
    end
end

end
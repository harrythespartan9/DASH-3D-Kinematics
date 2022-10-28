% This function returns an scalar error value based on the required
% amplification factor for each leg DOF, the error gradient, and hessian to
% support the fmincon approach to designing DASH kinematics. g_a and g_b
% are the amplifications needed in the swing and lift directions of the
% leg, and p contains details on which approximation of the amplification
% error to use.
function f = DASH_err_gamma(x, g_a, g_b, p) % , g, H add this to the output
                                              % when hessian is set

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

% Define the error here: 'l' intact
switch p.io

    case 0 % simpler error formulation based on IN-plane output components

        switch p.c
    
            case 0 % Doubly coupled mounting case
        
                f = (g_a^2 - (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))^2/(4*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2))^2 + ((2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))^2/(4*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2) - g_b^2)^2;
        
            case 1 % Swing lever, Lift coupled
                
                f = ((2*x(1)*cos(x(7)) + x(2)*cos(x(8)) - 2*x(4)*cos(x(10)) - 2*x(5)*cos(x(11)) - 2*x(6)*cos(x(12)) + 2*x(13)*cos(x(8)))^2/(4*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2) - g_a^2)^2 + (g_b^2 - (x(1)*sin(x(7) - pi/(12*x(14))) - x(1)*sin(x(7) + pi/(12*x(14))) + (x(2)*sin(x(8) - pi/(12*x(14))))/2 - (x(2)*sin(x(8) + pi/(12*x(14))))/2 + x(13)*sin(x(8) - pi/(12*x(14))) - x(13)*sin(x(8) + pi/(12*x(14))) + 2*sin(pi/(12*x(14)))*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))))^2/(4*sin(pi/(12*x(14)))^2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2))^2;
        
            case 2 % Swing coupled, Lift lever
        
                f = ((2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(11)))^2/(4*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2) - g_b^2)^2 + (g_a^2 - (x(4)*sin(x(10) - pi/(12*x(14))) - x(4)*sin(x(10) + pi/(12*x(14))) + (x(5)*sin(x(11) - pi/(12*x(14))))/2 - (x(5)*sin(x(11) + pi/(12*x(14))))/2 + x(13)*sin(x(11) - pi/(12*x(14))) - x(13)*sin(x(11) + pi/(12*x(14))) + 2*sin(pi/(12*x(14)))*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12))))^2/(4*sin(pi/(12*x(14)))^2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2))^2;
        
        end

    case 1 % complex error formulation based on OUT-of-plane components

        switch p.c
    
            case 0 % Doubly coupled mounting case
        
                f = (g_a^2 + ((2*cos(pi/(6*x(14))) - 2)*(4*x(6)^2*cos(x(12))^2 - 4*x(4)^2*cos(x(10))^2 + 4*x(13)^2*cos(x(8))^2 + 4*x(4)^2 + x(5)^2 + 2*x(5)*x(6)*cos(x(12) - x(11)) + 4*x(5)*x(13)*cos(x(8)) + 2*x(5)*x(6)*cos(x(11))*cos(x(12)) + 4*x(4)*x(5)*sin(x(10))*sin(x(11)) - 2*x(5)*x(6)*sin(x(11))*sin(x(12)) - 8*x(5)*x(13)*cos(x(8))*cos(x(11))^2 + 8*x(4)*x(13)*cos(x(8))*sin(x(10))*sin(x(11)) - 8*x(6)*x(13)*cos(x(8))*cos(x(11))*cos(x(12))))/(16*sin(pi/(12*x(14)))^2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2))^2 + (g_b^2 - (sin(pi/(12*x(14)))^2*(2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))^2 + sin(pi/(12*x(14)))^2*(2*x(1)*sin(x(7)) + x(2)*sin(x(8)) + 2*x(13)*sin(x(8)))^2)/(4*sin(pi/(12*x(14)))^2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2))^2;
        
            case 1 % Swing lever, Lift coupled
                
                f = (g_b^2 - ((x(1)*cos(x(7) - pi/(12*x(14))) - x(1)*cos(x(7) + pi/(12*x(14))) + (x(2)*cos(x(8) - pi/(12*x(14))))/2 - (x(2)*cos(x(8) + pi/(12*x(14))))/2 + x(13)*cos(x(8) - pi/(12*x(14))) - x(13)*cos(x(8) + pi/(12*x(14))))^2 + (x(1)*sin(x(7) - pi/(12*x(14))) - x(1)*sin(x(7) + pi/(12*x(14))) + (x(2)*sin(x(8) - pi/(12*x(14))))/2 - (x(2)*sin(x(8) + pi/(12*x(14))))/2 + x(13)*sin(x(8) - pi/(12*x(14))) - x(13)*sin(x(8) + pi/(12*x(14))) + 2*sin(pi/(12*x(14)))*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))))^2)/(4*sin(pi/(12*x(14)))^2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2))^2 + (g_a^2 + ((2*cos(pi/(6*x(14))) - 2)*(2*x(1)*cos(x(7)) + x(2)*cos(x(8)) - 2*x(4)*cos(x(10)) - 2*x(5)*cos(x(11)) - 2*x(6)*cos(x(12)) + 2*x(13)*cos(x(8)))^2)/(16*sin(pi/(12*x(14)))^2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2))^2;
        
            case 2 % Swing coupled, Lift lever
        
                f = (g_a^2 - ((x(4)*cos(x(10) - pi/(12*x(14))) - x(4)*cos(x(10) + pi/(12*x(14))) + (x(5)*cos(x(11) - pi/(12*x(14))))/2 - (x(5)*cos(x(11) + pi/(12*x(14))))/2 + x(13)*cos(x(11) - pi/(12*x(14))) - x(13)*cos(x(11) + pi/(12*x(14))))^2 + (x(4)*sin(x(10) - pi/(12*x(14))) - x(4)*sin(x(10) + pi/(12*x(14))) + (x(5)*sin(x(11) - pi/(12*x(14))))/2 - (x(5)*sin(x(11) + pi/(12*x(14))))/2 + x(13)*sin(x(11) - pi/(12*x(14))) - x(13)*sin(x(11) + pi/(12*x(14))) + 2*sin(pi/(12*x(14)))*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12))))^2)/(4*sin(pi/(12*x(14)))^2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2))^2 + (g_b^2 + ((2*cos(pi/(6*x(14))) - 2)*(2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(11)))^2)/(16*sin(pi/(12*x(14)))^2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2))^2;
        
        end

end

% % % % Define the error here: LEGACY CODE, NOT SUPPORTED
% % % f = (-(x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
% % %     (2*(x(5)*cos(x(11)) + x(4)*cos(x(10)) + x(6)*cos(x(12)))) - g_a)^2 + ... % swing squared error
% % %     ...
% % %     (-(2*x(2)*cos(x(8)) - x(5)*cos(x(11)) + 2*x(1)*cos(x(7)) - 2*x(4)*cos(x(10)) + ...
% % %     2*x(3)*cos(x(9)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
% % %     (2*(x(2)*cos(x(8)) + x(1)*cos(x(7)) + x(3)*cos(x(9)))) - g_b)^2; % lift squared error


% % Gradient and Hessian below: LEGACY CODE, NOT SUPPORTED % can still be
% % used as an example to setup the gradient and hessian if needed.
% if nargout > 1 % gradient required
%     
%     g = [];
% 
% %     if nargout > 2 % hessian required
% %         
% %         H = [];
% % 
% %     end
% end

end
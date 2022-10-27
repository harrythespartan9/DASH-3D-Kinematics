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

% % % % Define the error here: LEGACY CODE, NOT SUPPORTED
% % % f = (-(x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
% % %     (2*(x(5)*cos(x(11)) + x(4)*cos(x(10)) + x(6)*cos(x(12)))) - g_a)^2 + ... % swing squared error
% % %     ...
% % %     (-(2*x(2)*cos(x(8)) - x(5)*cos(x(11)) + 2*x(1)*cos(x(7)) - 2*x(4)*cos(x(10)) + ...
% % %     2*x(3)*cos(x(9)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
% % %     (2*(x(2)*cos(x(8)) + x(1)*cos(x(7)) + x(3)*cos(x(9)))) - g_b)^2; % lift squared error

% Define the error here: 'l' intact
switch p.io

    case 0 % simpler error formulation based on IN-plane output components

        switch p.c
    
            case 0 % Doubly coupled mounting case
        
                f = ((l_5*cos(theta_4 + theta_5) + 2*l_6*cos(theta_4 + theta_5 + theta_6) - 2*L*cos(theta_1 + theta_2)*cos(theta_4 + theta_5))^2/(4*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6))^2) - gamma_alpha^2)^2 + ((2*l_2*cos(theta_1 + theta_2) - l_5*cos(theta_4 + theta_5) + 2*l_1*cos(theta_1) - 2*l_4*cos(theta_4) + 2*l_3*cos(theta_1 + theta_2 + theta_3) - 2*L*cos(theta_1 + theta_2)*cos(theta_4 + theta_5))^2/(4*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3))^2) - gamma_beta^2)^2;
        
            case 1 % Swing lever, Lift coupled
                
                f = ((2*L*cos(theta_1 + theta_2) + l_2*cos(theta_1 + theta_2) - 2*l_5*cos(theta_4 + theta_5) + 2*l_1*cos(theta_1) - 2*l_4*cos(theta_4) - 2*l_6*cos(theta_4 + theta_5 + theta_6))^2/(4*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6))^2) - gamma_alpha^2)^2 + (gamma_beta^2 - ((l_2*sin(theta_1 + theta_2 - sym(pi)/(12*l)))/2 - (l_2*sin(theta_1 + theta_2 + sym(pi)/(12*l)))/2 + l_1*sin(theta_1 - sym(pi)/(12*l)) - l_1*sin(theta_1 + sym(pi)/(12*l)) + 2*sin(sym(pi)/(12*l))*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3)) + L*sin(theta_1 + theta_2 - sym(pi)/(12*l)) - L*sin(theta_1 + theta_2 + sym(pi)/(12*l)))^2/(4*sin(sym(pi)/(12*l))^2*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3))^2))^2;
        
            case 2 % Swing coupled, Lift lever
        
                f = ((2*L*cos(theta_4 + theta_5) - 2*l_2*cos(theta_1 + theta_2) + l_5*cos(theta_4 + theta_5) - 2*l_1*cos(theta_1) + 2*l_4*cos(theta_4) - 2*l_3*cos(theta_1 + theta_2 + theta_3))^2/(4*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3))^2) - gamma_beta^2)^2 + (gamma_alpha^2 - ((l_5*sin(theta_4 + theta_5 - sym(pi)/(12*l)))/2 - (l_5*sin(theta_4 + theta_5 + sym(pi)/(12*l)))/2 + l_4*sin(theta_4 - sym(pi)/(12*l)) - l_4*sin(theta_4 + sym(pi)/(12*l)) + 2*sin(sym(pi)/(12*l))*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6)) + L*sin(theta_4 + theta_5 - sym(pi)/(12*l)) - L*sin(theta_4 + theta_5 + sym(pi)/(12*l)))^2/(4*sin(sym(pi)/(12*l))^2*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6))^2))^2;
        
        end

    case 1 % complex error formulation based on OUT-of-plane components

        switch p.c
    
            case 0 % Doubly coupled mounting case
        
                f = (gamma_alpha^2 + ((2*cos(sym(pi)/(6*l)) - 2)*(l_5^2*cos(theta_4 + theta_5)^2 + l_5^2*sin(theta_4 + theta_5)^2 + 4*l_4^2*sin(theta_4)^2 + 4*l_6^2*cos(theta_4 + theta_5 + theta_6)^2 + 2*l_5*l_6*cos(theta_6) + 2*l_5*l_6*cos(2*theta_4 + 2*theta_5 + theta_6) + 4*L^2*cos(theta_1 + theta_2)^2*cos(theta_4 + theta_5)^2 + 4*L^2*cos(theta_1 + theta_2)^2*sin(theta_4 + theta_5)^2 - 4*L*l_5*cos(theta_1 + theta_2)*cos(theta_4 + theta_5)^2 + 4*l_4*l_5*sin(theta_4 + theta_5)*sin(theta_4) + 4*L*l_5*cos(theta_1 + theta_2)*sin(theta_4 + theta_5)^2 + 8*L*l_4*cos(theta_1 + theta_2)*sin(theta_4 + theta_5)*sin(theta_4) - 8*L*l_6*cos(theta_4 + theta_5 + theta_6)*cos(theta_1 + theta_2)*cos(theta_4 + theta_5)))/(16*sin(sym(pi)/(12*l))^2*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6))^2))^2 + (gamma_beta^2 - (sin(sym(pi)/(12*l))^2*(2*L*sin(theta_1 + theta_2) + l_2*sin(theta_1 + theta_2) + 2*l_1*sin(theta_1))^2 + sin(sym(pi)/(12*l))^2*(2*l_2*cos(theta_1 + theta_2) - l_5*cos(theta_4 + theta_5) + 2*l_1*cos(theta_1) - 2*l_4*cos(theta_4) + 2*l_3*cos(theta_1 + theta_2 + theta_3) - 2*L*cos(theta_1 + theta_2)*cos(theta_4 + theta_5))^2)/(4*sin(sym(pi)/(12*l))^2*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3))^2))^2;
        
            case 1 % Swing lever, Lift coupled
                
                f = (gamma_alpha^2 + ((2*cos(sym(pi)/(6*l)) - 2)*(2*L*cos(theta_1 + theta_2) + l_2*cos(theta_1 + theta_2) - 2*l_5*cos(theta_4 + theta_5) + 2*l_1*cos(theta_1) - 2*l_4*cos(theta_4) - 2*l_6*cos(theta_4 + theta_5 + theta_6))^2)/(16*sin(sym(pi)/(12*l))^2*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6))^2))^2 + (gamma_beta^2 - ((l_1*cos(theta_1 - sym(pi)/(12*l)) - l_1*cos(theta_1 + sym(pi)/(12*l)) + L*cos(theta_1 + theta_2 - sym(pi)/(12*l)) - L*cos(theta_1 + theta_2 + sym(pi)/(12*l)) + (l_2*cos(theta_1 + theta_2 - sym(pi)/(12*l)))/2 - (l_2*cos(theta_1 + theta_2 + sym(pi)/(12*l)))/2)^2 + ((l_2*sin(theta_1 + theta_2 - sym(pi)/(12*l)))/2 - (l_2*sin(theta_1 + theta_2 + sym(pi)/(12*l)))/2 + l_1*sin(theta_1 - sym(pi)/(12*l)) - l_1*sin(theta_1 + sym(pi)/(12*l)) + 2*sin(sym(pi)/(12*l))*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3)) + L*sin(theta_1 + theta_2 - sym(pi)/(12*l)) - L*sin(theta_1 + theta_2 + sym(pi)/(12*l)))^2)/(4*sin(sym(pi)/(12*l))^2*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3))^2))^2;
        
            case 2 % Swing coupled, Lift lever
        
                f = (gamma_beta^2 + ((2*cos(sym(pi)/(6*l)) - 2)*(2*L*cos(theta_4 + theta_5) - 2*l_2*cos(theta_1 + theta_2) + l_5*cos(theta_4 + theta_5) - 2*l_1*cos(theta_1) + 2*l_4*cos(theta_4) - 2*l_3*cos(theta_1 + theta_2 + theta_3))^2)/(16*sin(sym(pi)/(12*l))^2*(l_2*cos(theta_1 + theta_2) + l_1*cos(theta_1) + l_3*cos(theta_1 + theta_2 + theta_3))^2))^2 + (gamma_alpha^2 - ((l_4*cos(theta_4 - sym(pi)/(12*l)) - l_4*cos(theta_4 + sym(pi)/(12*l)) + L*cos(theta_4 + theta_5 - sym(pi)/(12*l)) - L*cos(theta_4 + theta_5 + sym(pi)/(12*l)) + (l_5*cos(theta_4 + theta_5 - sym(pi)/(12*l)))/2 - (l_5*cos(theta_4 + theta_5 + sym(pi)/(12*l)))/2)^2 + ((l_5*sin(theta_4 + theta_5 - sym(pi)/(12*l)))/2 - (l_5*sin(theta_4 + theta_5 + sym(pi)/(12*l)))/2 + l_4*sin(theta_4 - sym(pi)/(12*l)) - l_4*sin(theta_4 + sym(pi)/(12*l)) + 2*sin(sym(pi)/(12*l))*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6)) + L*sin(theta_4 + theta_5 - sym(pi)/(12*l)) - L*sin(theta_4 + theta_5 + sym(pi)/(12*l)))^2)/(4*sin(sym(pi)/(12*l))^2*(l_5*cos(theta_4 + theta_5) + l_4*cos(theta_4) + l_6*cos(theta_4 + theta_5 + theta_6))^2))^2;
        
        end

end



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
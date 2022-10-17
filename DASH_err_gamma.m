% This function returns an scalar error value based on the required
% amplification factor for each leg DOF, the error gradient, and hessian to
% support the fmincon approach to designing DASH kinematics.
function [f] = DASH_err_gamma(x, g_a, g_b) % , g, H add this to the output
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

% % % % Define the error here: OLD FORMAT
% % % f = (-(x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
% % %     (2*(x(5)*cos(x(11)) + x(4)*cos(x(10)) + x(6)*cos(x(12)))) - g_a)^2 + ... % swing squared error
% % %     ...
% % %     (-(2*x(2)*cos(x(8)) - x(5)*cos(x(11)) + 2*x(1)*cos(x(7)) - 2*x(4)*cos(x(10)) + ...
% % %     2*x(3)*cos(x(9)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
% % %     (2*(x(2)*cos(x(8)) + x(1)*cos(x(7)) + x(3)*cos(x(9)))) - g_b)^2; % lift squared error

% Define the error here: 'l' intact
f = (g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/...
    (2*x(14)))^2 + ...
    (g_b + ...
    (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - ...
    x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*x(14)))^2;

% % Gradient and Hessian below:
% if nargout > 1 % gradient required
%     
%     g = [(cos(x(7))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2; 
%         (cos(x(8))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2; 
%         (cos(x(9))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2; 
%         - (2*cos(x(10))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))) - (cos(x(10))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))))*(x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2; 
%         (cos(x(11))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))))*(x(4)*cos(x(10)) - x(6)*cos(x(12)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2 - (cos(x(11))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))); 
%         (cos(x(12))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2; 
%         -(x(1)*sin(x(7))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2; 
%         (2*x(13)*cos(x(11))*sin(x(8))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12))))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12))) - (sin(x(8))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(2*x(2)*x(4)*cos(x(10)) + x(2)*x(5)*cos(x(11)) - 2*x(1)*x(13)*cos(x(7))*cos(x(11)) - 2*x(3)*x(13)*cos(x(9))*cos(x(11))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2; 
%         -(x(3)*sin(x(9))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))^2; 
%         (2*x(4)*sin(x(10))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))) + (x(4)*sin(x(10))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))))*(x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2; 
%         (sin(x(11))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))))*(x(5) + 2*x(13)*cos(x(8))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))) + (sin(x(11))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))))*(x(5)*x(6)*cos(x(12)) - x(4)*x(5)*cos(x(10)) + 2*x(4)*x(13)*cos(x(8))*cos(x(10)) + 2*x(6)*x(13)*cos(x(8))*cos(x(12))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2; 
%         -(x(6)*sin(x(12))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))))*(2*x(4)*cos(x(10)) + x(5)*cos(x(11)) + 2*x(13)*cos(x(8))*cos(x(11))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12)))^2; 
%         - (2*cos(x(8))*cos(x(11))*(g_a + (x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12))))))/(x(4)*cos(x(10)) + x(5)*cos(x(11)) + x(6)*cos(x(12))) - (2*cos(x(8))*cos(x(11))*(g_b + (2*x(1)*cos(x(7)) + 2*x(2)*cos(x(8)) + 2*x(3)*cos(x(9)) - 2*x(4)*cos(x(10)) - x(5)*cos(x(11)) - 2*x(13)*cos(x(8))*cos(x(11)))/(2*(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9))))))/(x(1)*cos(x(7)) + x(2)*cos(x(8)) + x(3)*cos(x(9)))];
% 
% %     if nargout > 2 % hessian required
% %         
% %         H = [];
% % 
% %     end
% end

end
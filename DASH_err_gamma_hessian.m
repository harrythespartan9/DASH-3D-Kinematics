% This function returns an scalar error value based on the required
% amplification factor for each leg DOF
function f = DASH_err_gamma_hessian(x, g_a, g_b)

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

% Define the error here:
f = (-(x(5)*cos(x(11)) + 2*x(6)*cos(x(12)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
    (2*(x(5)*cos(x(11)) + x(4)*cos(x(10)) + x(6)*cos(x(12)))) - g_a)^2 + ... % swing squared error
    ...
    (-(2*x(2)*cos(x(8)) - x(5)*cos(x(11)) + 2*x(1)*cos(x(7)) - 2*x(4)*cos(x(10)) + ...
    2*x(3)*cos(x(9)) - 2*x(13)*cos(x(8))*cos(x(11)))/ ...
    (2*(x(2)*cos(x(8)) + x(1)*cos(x(7)) + x(3)*cos(x(9)))) - g_b)^2; % lift squared error

end
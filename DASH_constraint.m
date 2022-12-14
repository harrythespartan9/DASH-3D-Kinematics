% This function returns an scalar error value based on the required
% amplification factor for each leg DOF
function [c, ceq] = DASH_constraint(x, ~, ~) % , gradc, gradceq

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

% We have no inequality constraints, so set that as empty.
c = [];

% % % % Define the equality constraints here: OLD FORMAT
% % % ceq(1) = x(2)*sin(x(8)) + x(1)*sin(x(7)) + x(3)*sin(x(9));
% % % ceq(2) = -x(5)*sin(x(11)) -x(4)*sin(x(10)) -x(6)*sin(x(12));
% % % ceq(3) = x(2)*cos(x(8)) - x(5)*cos(x(11)) + ...
% % %     x(1)*cos(x(7)) - x(4)*cos(x(10)) + ...
% % %     x(3)*cos(x(9)) - x(6)*cos(x(12));

% Define the equality constraints here:
ceq(1) = x(1)*sin(x(7)) + x(2)*sin(x(8)) + x(3)*sin(x(9));
ceq(2) = x(4)*sin(x(10)) + x(5)*sin(x(11)) + x(6)*sin(x(12));
ceq(3) = x(14) - x(1)*cos(x(7)) - x(2)*cos(x(8)) - x(3)*cos(x(9));
ceq(4) = x(14) - x(4)*cos(x(10)) - x(5)*cos(x(11)) - x(6)*cos(x(12));

% The constraints above satisfy the coupled kinematic length and lateral
% length violations. Refer to line numbers 110, 111, 112, and 113 in 
% "DASH_leg_kinematics.mlx" for more details.

% % If constraint gradients are requested, get the outputs we need.
% if nargout > 2
%     gradc = [];
%     gradceq = [[sin(x(7)); sin(x(8)); sin(x(9)); 0; 0; 0; x(1)*cos(x(7)); x(2)*cos(x(8)); x(3)*cos(x(9)); 0; 0; 0; 0],...
%         [0; 0; 0; sin(x(10)); sin(x(11)); sin(x(12)); 0; 0; 0; x(4)*cos(x(10)); x(5)*cos(x(11)); x(6)*cos(x(12)); 0],...
%         [cos(x(7)); cos(x(8)); cos(x(9)); -cos(x(10)); -cos(x(11)); -cos(x(12)); -x(1)*sin(x(7)); -x(2)*sin(x(8)); -x(3)*sin(x(9)); x(4)*sin(x(10)); x(5)*sin(x(11)); x(6)*sin(x(12)); 0]];
% end

end
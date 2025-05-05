function joint_values = InverseKinematicsRPR(l1, l2, x, y, alpha, pos_neg)
    % Function that performs the inverse kinematic of an RPR robot
    % Parameters:
    % - l2: length of the second link
    % - x: x coordinate of end effector position
    % - y: y coordinate of end effector position
    % - alpha: angle of the end effector (most of the time given by q1 + q3)
    % - pos_neg: a string that indicates what solution you want "pos" or "neg"

    % Calculate the joint values
    if pos_neg == "pos"
        q2 = sqrt((x - l2 * cos(alpha))^2 + (y - l2 * sin(alpha))^2 - l1^2);
    else
        q2 = -sqrt((x - l2 * cos(alpha))^2 + (y - l2 * sin(alpha))^2 - l1^2);
    end

    q1 = atan2(l1*(y - l2 * sin(alpha)) - q2*(x - l2*cos(alpha)), l1*(x - l2 * cos(alpha)) + q2*(y - l2*sin(alpha)));
    q3 = alpha - q1;

    disp("Joint values:")
    disp("- q1:");
    disp(q1);
    disp("- q2:");
    disp(q2);
    disp("- q3:");
    disp(q3);

    joint_values = [q1, q2, q3];
end

function dtheta_val = dtheta(theta1, theta2)
% in periodic [0, pi].
% Same size, or matrix / 1, or 1 / matrix.

if (length(theta1(:)) ~= 1) && (length(theta2(:)) == 1)
    theta2 = ones(size(theta1)) * theta2;
elseif (length(theta1(:)) == 1) && (length(theta2(:)) ~= 1)
    theta1 = ones(size(theta2)) * theta1;
elseif (length(theta1(:)) ~= 1) && (length(theta2(:)) ~= 1) && (~isequal(size(theta1), size(theta2)))
    error('theta1 and theta2 are not in the same size.');
end

% inner product
dtheta_abs = acos(cos(2 * theta1) .* cos(2 * theta2) + sin(2 * theta1) .* sin(2 * theta2)) / 2;
% cross product
dtheta_sgn = sign(cos(2 * theta2) .* sin(2 * theta1) - cos(2 * theta1) .* sin(2 * theta2));
dtheta_val = dtheta_abs .* dtheta_sgn;



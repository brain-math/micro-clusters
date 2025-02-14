function c_periodic_fft = conv2_periodic_by_fft2(x, kernel)

% https://www.mathworks.com/help/signal/ug/linear-and-circular-convolution.html
% If x and y have the same size, then ifft of the product of fft of x and y is the result of circular convolution.
% But there are two problems -- 1) sometimes we only want x but not y to be periodic, eg. x is 2D neuron population, y is kernel;
    % 2) x and y could be in different size...
% Solution: manually add periodic padding to x, then do linear conv.
% To do linear conv, zero pad both x and y then fft is equiv. to conv with 'full' (default). Then cut to the size of x.
%
% BTW, if x and kernel have the same size, then they're commutative.

m = size(kernel);
if ~isequal(mod(m, 2), [1, 1]), error('Kernel size should be odd.'); end

% Periodic padding x, so circular BC of x goes back to a linear conv. problem.
na = (m - 1) / 2;
xpad = x;
xpad = cat(1, xpad(end - (na(1) - 1): end, :), xpad, xpad(1: na(1), :));
xpad = cat(2, xpad(:, end - (na(2) - 1): end), xpad, xpad(:, 1: na(2)));
%
% zero pad both x and y then fft is equiv. to conv with 'full' (default). 
L = size(xpad) + m - 1;
xpad2 = zeros(L); xpad2(1: size(xpad, 1), 1: size(xpad, 2)) = xpad;
ypad = zeros(L); ypad(1: m(1), 1: m(2)) = kernel;
c_periodic_fft = ifft2(fft2(xpad2) .* fft2(ypad));
%
% Cut from 'full' to the size of input x.
c_periodic_fft = c_periodic_fft(L(1) - size(xpad, 1) + 1: size(xpad, 1),...
    L(2) - size(xpad, 2) + 1: size(xpad, 2));


function [x, n_iter, error] = conjGradFunc(A, b, x0, m, er_max)
    % Conjugate Gradient algorithm
    w = zeros(size(A, 1), m+1);
    p = zeros(size(A, 1), m+1);
    alpha = zeros(1, m+1);
    r = zeros(size(A, 1), m+1);
    r(:, 1) = -A*x0 + b;
    w(:, 1) = r(:, 1);
    for j = 2:m+1
        p(:, j-1) = A*w(:, j-1);
        alpha(j-1) = r(:, j-1).' * r(:, j-1) / (p(:, j-1).' * r(:, j-1));
        x0 = x0 + alpha(j-1)*w(:, j-1);
        r(:, j) = r(:, j-1) - alpha(j-1)*p(:, j-1);
        error = norm(r(:, j));
        beta = r(:, j).' * r(:, j) / (r(:, j-1).' * r(:, j-1));
        w(:, j) = r(:, j) + beta*w(:, j-1);
        if error < er_max
            break;
        end
    end
    x = x0;
    n_iter = j-1;
end

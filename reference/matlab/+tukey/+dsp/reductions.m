function g = reductions()
    % REDUCTIONS Reductions: sum and mean, differences and gradients,
    % running sums and trapezoidal integrals, and percentiles.
    categories = {"sum", "mean", "diff", "gradient", "cumsum", ...
                  "trapezoid", "cumulative_trapezoid", "percentile"};
    runs = {@sumVector, @meanVector, @diffVector, @gradientVector, ...
            @cumsumVector, @trapzVector, @cumtrapzVector, @prctileVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = sumVector(params)
    output = struct("sum", encode.vec(sum(params.x)));
    description = "sum(x)";
end

function [output, description] = meanVector(params)
    output = struct("mean", encode.vec(mean(params.x)));
    description = "mean(x)";
end

function [output, description] = diffVector(params)
    output = struct("d", encode.vec(diff(params.x)));
    description = "diff(x)";
end

function [output, description] = gradientVector(params)
    output = struct("g", encode.vec(gradient(params.x, params.dx)));
    description = "gradient(x, dx)";
end

function [output, description] = cumsumVector(params)
    output = struct("s", encode.vec(cumsum(params.x)));
    description = "cumsum(x)";
end

function [output, description] = trapzVector(params)
    % A scalar first argument is the sample spacing. The dimension is given
    % explicitly, as with two scalar arguments (a one-sample y) trapz reads
    % trapz(Y, dim) instead of trapz(X, Y).
    output = struct("integral", ...
                    encode.vec(trapz(params.dx, params.y(:), 1)));
    description = "trapz(dx, y, 1)";
end

function [output, description] = cumtrapzVector(params)
    % As trapz, the dimension is explicit for a one-sample y.
    output = struct("integral", ...
                    encode.vec(cumtrapz(params.dx, params.y(:), 1)));
    description = "cumtrapz(dx, y, 1)";
end

function [output, description] = prctileVector(params)
    % prctile's default method is checked against NumPy's "hazen". It has
    % no equivalent of NumPy's default, "linear", so those cases are
    % skipped. The NaN case is kept, to show how prctile treats NaN.
    output = [];
    description = "";
    if params.method ~= "hazen"
        return
    end
    output = struct("values", encode.vec(prctile(params.x, params.p)));
    description = "prctile(x, p)";
end

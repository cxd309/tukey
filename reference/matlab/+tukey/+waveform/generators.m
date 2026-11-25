function g = generators()
    % GENERATORS The tukey/waveform categories MATLAB can reproduce, and how
    % each run function takes a SciPy vector's params and returns MATLAB's
    % output (encoded with the encode namespace) and a description,
    % or an empty output when MATLAB has no equivalent for that case
    %
    % chirp, from MATLAB's chirp.m in R2025a:
    % - linear is SciPy's formula
    % - quadratic with its default shape, "concave" for an up sweep and
    %   "convex" for a down sweep, has its vertex at t = 0, as SciPy's
    %   vertex_zero=True
    % - logarithmic is SciPy's formula
    %
    % square, from MATLAB's square.m in R2025a:
    % - duty is a percentage, so SciPy's duty is multiplied by 100
    % - mod(t, 2*pi) rounds near whole multiples of 2*pi, where SciPy's % is
    %   exact, so the two may differ at edges
    %
    % sawtooth, from MATLAB's sawtooth.m in R2025a:
    % - width is a fraction of the period, as SciPy's
    % - t is reduced with rem, which keeps t's sign, then scaled to a fraction
    %   of a cycle, so a tiny negative t at width 1 gives 1, where SciPy's %
    %   rounds it up to 2*pi and gives NaN
    %
    % not reproduced:
    % - hyperbolic: no MATLAB method
    % - quadratic with its vertex at t1: MATLAB's other shape evaluates the
    %   sweep at fliplr(-t) with f0 and f1 swapped, so its output depends on
    %   t's span, order and orientation. On a row vector from 0 to t1 it is
    %   SciPy's sweep less a constant phase, and on a column it's reversed
    % - logarithmic with f0 = f1: MATLAB divides by log(f1/f0) = 0
    categories = {"waveform/chirp", "waveform/square", "waveform/sawtooth"};
    runs = {@chirpVector, @squareVector, @sawtoothVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = chirpVector(params)
    output = [];
    description = "";
    if params.method == "hyperbolic" || ~params.vertex_zero || ...
       (params.method == "logarithmic" && params.f0 == params.f1)
        return
    end
    x = chirp(params.t, params.f0, params.t1, params.f1, params.method, ...
              params.phi);
    output = struct("x", encode.vec(x));
    description = sprintf("chirp(t, %.17g, %.17g, %.17g, ""%s"", %.17g)", ...
                          params.f0, params.t1, params.f1, params.method, ...
                          params.phi);
end

function [output, description] = squareVector(params)
    % square takes its duty cycle as a percentage, where SciPy's is a fraction
    duty = 100 * params.duty;
    x = square(params.t, duty);
    output = struct("x", encode.vec(x));
    description = sprintf("square(t, %.17g)", duty);
end

function [output, description] = sawtoothVector(params)
    x = sawtooth(params.t, params.duty);
    output = struct("x", encode.vec(x));
    description = sprintf("sawtooth(t, %.17g)", params.duty);
end

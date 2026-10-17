function g = generators()
% GENERATORS The tukey/window categories MATLAB can reproduce, and how
% each run function takes a SciPy vector's params and returns MATLAB's
% output (encoded with the encode namespace) and a description,
% or an empty output when MATLAB has no equivalent for that case
%
% MATLAB's names and parameters differ from SciPy's:
% - boxcar is rectwin
% - nuttall, flattop, barthann, bohman and parzen add "win"
% - tukey is tukeywin (alpha as its r)
% - gaussian is gausswin, whose alpha is (n-1)/(2*std)
% - taylor is taylorwin, whose sll is negative dB
% - chebwin's r is SciPy's at
%
% windows taking "symmetric" or "periodic" are reproduced both ways
% the others only have a symmetric form, so their periodic cases are skipped
%
% hanning is MATLAB's own: SciPy's vectors build it from MATLAB's documented
% definition, which these check
%
% not reproduced (no MATLAB function): 
% - cosine
% - lanczos
% - exponential
% - general_cosine
% - general_hamming
% - general_gaussian
% - kaiser_bessel_derived
%
% taylorwin doesn't normalise, so only SciPy's norm=False cases are reproduced,
% SciPy returns [1] for every window there where taylorwin evaluates its formula
g = struct( ...
    "category", {"window/boxcar", "window/triang", "window/bartlett", ...
                 "window/hann", "window/hamming", "window/blackman", ...
                 "window/blackmanharris", "window/nuttall", "window/flattop", ...
                 "window/hanning", ...
                 "window/barthann", "window/bohman", "window/parzen", ...
                 "window/tukey", "window/kaiser", "window/gaussian", ...
                 "window/taylor", "window/chebwin"}, ...
    "run", {symmetricOnly(@rectwin), symmetricOnly(@triang), symmetricOnly(@bartlett), ...
            withSflag(@hann), withSflag(@hamming), withSflag(@blackman), ...
            withSflag(@blackmanharris), withSflag(@nuttallwin), withSflag(@flattopwin), ...
            withSflag(@hanning), ...
            symmetricOnly(@barthannwin), symmetricOnly(@bohmanwin), symmetricOnly(@parzenwin), ...
            symmetricOnly(@tukeywin, @(p) {p.alpha}), ...
            symmetricOnly(@kaiser, @(p) {p.beta}), ...
            symmetricOnly(@gausswin, @(p) {(p.n - 1) / (2 * p.std)}), ...
            onlyWhen(@(p) ~p.norm && p.n ~= 1, ...
                symmetricOnly(@taylorwin, @(p) {p.nbar, -p.sll})), ...
            symmetricOnly(@chebwin, @(p) {p.at})});
end

% -------------------------------------------------------------------------
% windows with a "symmetric" or "periodic" flag
% -------------------------------------------------------------------------

function run = withSflag(window)
% withSflag is the run function for a window called as window(n, sflag)
run = @(params) withSflagVector(window, params);
end

function [output, description] = withSflagVector(window, params)
sflag = "symmetric";
if ~params.sym
    sflag = "periodic";
end
w = window(params.n, sflag);
output = struct("w", encode.vec(w));
description = sprintf("%s(%d, ""%s"")", func2str(window), params.n, sflag);
end

% -------------------------------------------------------------------------
% windows with only a symmetric form
% -------------------------------------------------------------------------

function run = symmetricOnly(window, parameters)
% symmetricOnly is the run function for a window called as window(n, args{:}),
% where args = parameters(params) maps SciPy's parameters to MATLAB's,
% e.g. tukeywin's r from SciPy's alpha; windows without parameters get none
arguments
    window function_handle
    parameters function_handle = @(p) {}
end
run = @(params) symmetricOnlyVector(window, parameters, params);
end

function [output, description] = symmetricOnlyVector(window, parameters, params)
output = [];
description = "";
if ~params.sym
    return
end
args = parameters(params);
w = window(params.n, args{:});
output = struct("w", encode.vec(w));
description = sprintf("%s(%s)", func2str(window), ...
    strjoin(compose("%.17g", [params.n, args{:}]), ", "));
end

% -------------------------------------------------------------------------
% case filters
% -------------------------------------------------------------------------

function run = onlyWhen(predicate, inner)
% onlyWhen is inner's run function, skipping the cases where predicate(params) is false
run = @(params) onlyWhenVector(predicate, inner, params);
end

function [output, description] = onlyWhenVector(predicate, inner, params)
output = [];
description = "";
if predicate(params)
    [output, description] = inner(params);
end
end
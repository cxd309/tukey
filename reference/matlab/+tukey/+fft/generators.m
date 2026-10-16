function g = generators()
% GENERATORS The tukey/fft categories MATLAB can reproduce, and how
% each run function takes a SciPy vector's params and returns MATLAB's
% output (encoded with the encode namespace) and a description,
% or an empty output when MATLAB has no equivalent for that case
%
% MATLAB has no real-input transforms, so rfft is fft cropped to its
% first n/2+1 bins, and irfft is ifft with "symmetric"
%
% not reproduced: fftfreq (no MATLAB function), and the ortho and forward
% cases of fft_options (MATLAB has no normalisation option)
g = struct( ...
    "category", {"fft", "ifft", "rfft", "irfft", "fft_options", "fftshift"}, ...
    "run", {@fftVector, @ifftVector, @rfftVector, @irfftVector, ...
            @fftOptionVector, @fftshiftVector});
end

% -------------------------------------------------------------------------
% transforms
% -------------------------------------------------------------------------

function [output, description] = fftVector(params)
y = fft(toComplex(params.x));
output = struct("y", encode.complex(y));
description = "fft(x)";
end

function [output, description] = ifftVector(params)
y = ifft(toComplex(params.x));
output = struct("y", encode.complex(y));
description = "ifft(x)";
end

function [output, description] = rfftVector(params)
y = realForward(params.x, []);
output = struct("y", encode.complex(y));
description = "fft(x), first n/2+1 bins";
end

function [output, description] = irfftVector(params)
y = ifft(toComplex(params.x), params.n, "symmetric");
output = struct("y", encode.vec(y));
description = sprintf("ifft(x, %d, ""symmetric"")", params.n);
end

% -------------------------------------------------------------------------
% options and helpers
% -------------------------------------------------------------------------

function [output, description] = fftOptionVector(params)
% MATLAB has no norm option: only SciPy's default, "backward", has an equivalent
% an empty n leaves the length unchanged, in MATLAB as in SciPy
output = [];
description = "";
if params.norm ~= "backward"
    return
end
x = toComplex(params.x);
n = params.n;
switch params.transform
    case "fft"
        y = fft(x, n);
    case "ifft"
        y = ifft(x, n);
    case "rfft"
        y = realForward(real(x), n);
    case "irfft"
        y = ifft(x, n, "symmetric");
end
output = struct("y", encode.complex(y));
description = sprintf("%s(x, n=%s), as MATLAB", params.transform, mat2str(n));
end

function [output, description] = fftshiftVector(params)
output = struct("shift", encode.vec(fftshift(params.x)), ...
    "ishift", encode.vec(ifftshift(params.x)));
description = "fftshift(x) and ifftshift(x)";
end

function y = realForward(x, n)
% realForward is scipy.fft.rfft: fft(x, n) cropped to its first n/2+1 bins
if isempty(n)
    n = numel(x);
end
y = fft(x, n);
y = y(1:floor(n/2) + 1);
end

function x = toComplex(c)
% toComplex combines a decoded {"re": [...], "im": [...]} into complex values
x = complex(c.re, c.im);
end
function g = convolution()
    % CONVOLUTION Convolution and correlation: conv in its three shapes,
    % and xcorr.
    categories = {"convolve", "correlate"};
    runs = {@convVector, @xcorrVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = convVector(params)
    % 'same' is the central part the length of a, and 'valid' is empty
    % when v is longer than a, where SciPy swaps the inputs.
    a = params.a(:);
    v = params.v(:);
    output = struct("full", encode.vec(conv(a, v)), ...
                    "same", encode.vec(conv(a, v, "same")), ...
                    "valid", encode.vec(conv(a, v, "valid")));
    description = "conv(a, v), conv(a, v, ""same""), conv(a, v, ""valid"")";
end

function [output, description] = xcorrVector(params)
    % xcorr pads the shorter input with zeros to the longer's length and
    % gives every lag from -(N - 1) to N - 1. The lags outside SciPy's full
    % range, -(numel(v) - 1) to numel(a) - 1, sum only padding, so they're
    % left out, and the lags are written to check xcorr's direction.
    % xcorr has no 'same' or 'valid'. The maximum lag is given explicitly,
    % its default N - 1, as with a scalar v xcorr(a, v) reads v as the
    % maximum lag.
    a = params.a(:);
    v = params.v(:);
    [r, lags] = xcorr(a, v, max(numel(a), numel(v)) - 1);
    keep = lags >= -(numel(v) - 1) & lags <= numel(a) - 1;
    output = struct("full", encode.vec(r(keep)), ...
                    "full_lags", encode.vec(lags(keep)));
    description = "[r, lags] = xcorr(a, v, N - 1), at SciPy's full lags";
end

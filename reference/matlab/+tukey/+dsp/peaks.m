function g = peaks()
    % PEAKS Peak finding: findpeaks with each option, and the prominences
    % and widths of the peaks it finds.
    categories = {"find_peaks", "peak_prominences", "peak_widths"};
    runs = {@findPeaksVector, @peakProminencesVector, @peakWidthsVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = findPeaksVector(params)
    % SciPy's options, given as findpeaks' name-value pairs.
    % findpeaks gives 1-based locations, written 0-based to match SciPy.
    % It rejects signals of fewer than 3 samples, where SciPy finds no
    % peaks, and a MinPeakDistance of n - 1 or more, which SciPy accepts,
    % so those cases have no equivalent.
    output = [];
    description = "";
    n = numel(params.x);
    if n < 3 || (isfield(params, "distance") && params.distance >= n - 1)
        return
    end
    scipy = ["distance", "height", "threshold", "prominence", "width"];
    matlab = ["MinPeakDistance", "MinPeakHeight", "Threshold", ...
              "MinPeakProminence", "MinPeakWidth"];
    args = {};
    parts = "x";
    for i = 1:numel(scipy)
        if isfield(params, scipy(i))
            value = params.(scipy(i));
            args = [args, {matlab(i), value}]; %#ok<AGROW>
            parts(end + 1) = sprintf("""%s"", %g", ...
                                     matlab(i), value); %#ok<AGROW>
        end
    end
    [pks, locs] = findpeaks(params.x, args{:});
    output = struct("peaks", encode.vec(locs - 1), ...
                    "heights", encode.vec(pks));
    description = "findpeaks(" + strjoin(parts, ", ") + ")";
end

function [output, description] = peakProminencesVector(params)
    % findpeaks gives the prominence of each peak it finds, but not the
    % bases. Its peaks must be SciPy's for the comparison to mean anything.
    [~, locs, ~, p] = findpeaks(params.x);
    assert(isequal(locs(:) - 1, params.peaks(:)), ...
           "findpeaks found different peaks from SciPy");
    output = struct("prominences", encode.vec(p));
    description = "[~, ~, ~, p] = findpeaks(x)";
end

function [output, description] = peakWidthsVector(params)
    % findpeaks measures width at half prominence, SciPy's rel_height 0.5.
    % It has no equivalent at other heights.
    output = [];
    description = "";
    if params.rel_height ~= 0.5
        return
    end
    [~, locs, w] = findpeaks(params.x);
    assert(isequal(locs(:) - 1, params.peaks(:)), ...
           "findpeaks found different peaks from SciPy");
    output = struct("widths", encode.vec(w));
    description = "[~, ~, w] = findpeaks(x)";
end

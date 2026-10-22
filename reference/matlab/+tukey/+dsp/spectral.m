function g = spectral()
    % SPECTRAL Spectral analysis:
    % spectrogram's STFT, at FFT bins and at chosen frequencies, and its power,
    % pwelch and periodogram.
    % Every case passes its window as a vector and its FFT length explicitly,
    % as MATLAB's own defaults differ from SciPy's.
    % Time-frequency outputs are transposed to [segment][frequency],
    % tukey's layout.
    categories = {"stft", "stft_frequencies", "spectrogram", "welch", ...
                  "periodogram"};
    runs = {@stftVector, @stftFrequenciesVector, @spectrogramVector, ...
            @welchVector, @periodogramVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = stftVector(params)
    nfft = fftLength(params);
    if params.two_sided
        [s, f, t] = spectrogram(params.x, params.window, params.overlap, ...
                                nfft, params.fs, "twosided");
        description = ...
            "[S, F, T] = spectrogram(x, w, noverlap, nfft, fs, ""twosided"")";
    else
        [s, f, t] = spectrogram(params.x, params.window, params.overlap, ...
                                nfft, params.fs);
        description = "[S, F, T] = spectrogram(x, w, noverlap, nfft, fs)";
    end
    output = struct("frequencies", encode.vec(f), "times", encode.vec(t), ...
                    "values", complexRows(s.'));
end

function [output, description] = stftFrequenciesVector(params)
    % A vector in nfft's place evaluates the transform at those frequencies.
    [s, f, t] = spectrogram(params.x, params.window, params.overlap, ...
                            params.frequencies, params.fs);
    output = struct("frequencies", encode.vec(f), "times", encode.vec(t), ...
                    "values", complexRows(s.'));
    description = "[S, F, T] = spectrogram(x, w, noverlap, f, fs)";
end

function [output, description] = spectrogramVector(params)
    % The fourth output is the power, "psd" for SciPy's density and "power"
    % for its spectrum scaling.
    type = spectrumType(params.scaling);
    [~, f, t, p] = spectrogram(params.x, params.window, params.overlap, ...
                               fftLength(params), params.fs, type);
    output = struct("frequencies", encode.vec(f), "times", encode.vec(t), ...
                    "values", encode.rows(p.'));
    description = "[~, F, T, P] = spectrogram(x, w, noverlap, " + ...
                  "nfft, fs, """ + type + """)";
end

function [output, description] = welchVector(params)
    % pwelch averages with the mean only, and doesn't detrend each segment,
    % so those cases have no equivalent.
    output = [];
    description = "";
    if params.average ~= "mean" || params.detrend ~= "none"
        return
    end
    type = spectrumType(params.scaling);
    [p, f] = pwelch(params.x, params.window, params.overlap, ...
                    fftLength(params), params.fs, type);
    output = struct("frequencies", encode.vec(f), "power", encode.vec(p));
    description = "[pxx, f] = pwelch(x, w, noverlap, nfft, fs, """ + ...
                  type + """)";
end

function [output, description] = periodogramVector(params)
    type = spectrumType(params.scaling);
    [p, f] = periodogram(params.x, params.window, fftLength(params), ...
                         params.fs, type);
    output = struct("frequencies", encode.vec(f), "power", encode.vec(p));
    description = sprintf("[pxx, f] = periodogram(x, w, nfft, fs, ""%s"")", ...
                          type);
end

function nfft = fftLength(params)
    % The window's length when the SciPy case leaves nfft out (null), as
    % MATLAB would otherwise round it up to a power of two of at least 256.
    nfft = params.nfft;
    if isempty(nfft)
        nfft = numel(params.window);
    end
end

function type = spectrumType(scaling)
    % MATLAB's name for SciPy's scaling.
    if scaling == "spectrum"
        type = "power";
    else
        type = "psd";
    end
end

function s = complexRows(m)
    % complexRows encodes a complex matrix as {"re": rows, "im": rows}.
    s = "{""re"": " + encode.rows(real(m)) + ", ""im"": " + ...
        encode.rows(imag(m)) + "}";
end

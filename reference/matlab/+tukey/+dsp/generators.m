function g = generators()
    % GENERATORS The tukey/dsp categories MATLAB can reproduce, and how
    % each run function takes a SciPy vector's params and returns MATLAB's
    % output (encoded with the encode namespace) and a description,
    % or an empty output when MATLAB has no equivalent for that case.
    %
    % One file per family, as reference/scipy/dsp, each returning its own
    % categories and run functions.
    g = [tukey.dsp.design(), tukey.dsp.filtering(), tukey.dsp.response(), ...
         tukey.dsp.conditioning(), tukey.dsp.peaks(), tukey.dsp.envelopes()];
end

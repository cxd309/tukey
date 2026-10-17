function generate(scipyDir, outDir)
% GENERATE Write tukey's MATLAB reference vectors
% generate(scipyDir, outDir) runs the MATLAB equivalent of each SciPy
% reference vector in scipyDir and writes outDir/<category>/<name>.json
%
% Each file copies its SciPy vector's params verbatim, so both sources
% are checked on identical inputs, and replaces the output with MATLAB's
%
% Generators live in one namespace per tukey Go package (+tukey/+dsp, ...),
% and encode their outputs with the encode namespace
%
% Run with `just vectors-matlab` which clears outDir first

arguments
    scipyDir (1,1) string
    outDir (1,1) string
end

signal = ver("signal");
versions = struct("matlab", version, "signal", signal.Version);

% every generator, package by package, as in reference/scipy/main.py
generators = [tukey.dsp.generators(), tukey.fft.generators(), tukey.window.generators()];
for g = generators
    eachVector(scipyDir, outDir, versions, g.category, g.run);
end
end

% -------------------------------------------------------------------------
% reading and writing vectors
% -------------------------------------------------------------------------

function eachVector(scipyDir, outDir, versions, category, run)
% eachVector runs every SciPy vector in category through run and writes
%   run(params) returns the output as a struct of encoded JSON values
%   (see the encode namespace) and a description, an empty output skips the case
files = dir(fullfile(scipyDir, category, "*.json"));
n = 0;
for f = files'
    text = string(fileread(fullfile(f.folder, f.name)));
    v = jsondecode(text);
    [output, description] = run(v.params);
    if isempty(output)
        continue
    end
    writeVector(fullfile(outDir, category, f.name), description, ...
        versions, rawParams(text), output);
    n = n + 1;
end
fprintf("%s: %d\n", category, n);
end

function params = rawParams(text)
% rawParams is the SciPy vector's "params" exactly as written
% relies on reference/scipy/write.py's layout: indent 2, "output" last
params = extractBetween(text, """params"": ", ...
    "," + newline + "  ""output"": ");
assert(isscalar(params), "can't find params in a SciPy vector")
end

function writeVector(path, description, versions, params, output)
% writeVector writes one vector; output's values are already encoded JSON
entries = compose("""%s"": %s", string(fieldnames(output)), ...
    string(struct2cell(output)));
text = sprintf("{""description"": %s, ""versions"": %s, " + ...
    """params"": %s, ""output"": {%s}}", jsonencode(description), ...
    jsonencode(versions), params, strjoin(entries, ", "));

folder = fileparts(path);
if ~isfolder(folder)
    mkdir(folder);
end
writelines(text, path);
end

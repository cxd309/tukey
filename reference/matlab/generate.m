function generate(scipyDir, outDir)
%GENERATE Write tukey's MATLAB reference vectors
%   generate(scipyDir, outDir) runs the MATLAB equivalent of each SciPy
%   reference vector in scipyDir and writes outDir/<category>/<name>.json
%
%   Each file copies its SciPy vector's params verbatim, so both sources
%   are checked on identical inputs, and replaces the output with MATLAB's
%
%   Run with `just vectors-matlab` which clears outDir first

arguments
    scipyDir (1,1) string
    outDir (1,1) string
end

signal = ver("signal");
versions = struct("matlab", version, "signal", signal.Version);

eachVector(scipyDir, outDir, versions, "butter", @butterVector);
end

% -------------------------------------------------------------------------
% dsp
% -------------------------------------------------------------------------

function [output, description] = butterVector(params)
% butter's (b, a), and its zpk through zp2sos for the sos output
bandTypes = dictionary(...
    ["lowpass", "highpass", "bandpass", "bandstop"], ...
    ["low", "high", "bandpass", "stop"]);
btype = bandTypes(params.btype);
wn = params.wn';

[b, a] = butter(params.order, wn, btype);
[z, p, k] = butter(params.order, wn, btype);
sos = zp2sos(z, p, k);

output = struct("b", vec(b), "a", vec(a), "sos", rows(sos));
description = sprintf("butter(%d, %s, ""%s""); zp2sos of its zpk", ...
    params.order, mat2str(wn), btype);
end

% -------------------------------------------------------------------------
% reading and writing vectors
% -------------------------------------------------------------------------

function eachVector(scipyDir, outDir, versions, category, run)
% eachVector runs every SciPy vector in category through run and writes
%   run(params) returns the output as a struct of encoded JSON values
%   (see vec, rows, cvec) and a description, an empty output skips the case
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

function s = vec(x)
% vec encodes a real vector as a JSON array at full precision
% works on a single row
s = "[" + strjoin(compose("%.17g", x(:)'), ", ") + "]";
end

function s = rows(m)
% rows encodes a matrix as a JSON array of its rows
% works on a single row
r = strings(size(m, 1), 1);
for i = 1:size(m, 1)
    r(i) = vec(m(i, :));
end
s = "[" + strjoin(r, ", ") + "]";
end

function s = cvec(z)
% cvec encodes a complex vector as {"re": [...], "im": [...]}
s = "{""re"": " + vec(real(z)) + ", ""im"": " + vec(imag(z)) + "}";
end
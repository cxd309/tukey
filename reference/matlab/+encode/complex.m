function s = complex(z)
    % COMPLEX Encode a complex vector as {"re": [...], "im": [...]}
    re = encode.vec(real(z));
    im = encode.vec(imag(z));
    s = "{""re"": " + re + ", ""im"": " + im + "}";
end

testdata := justfile_directory() / "testdata"
scipy_reference_dir := "reference/scipy"
matlab := "/Applications/MATLAB_R2025a.app/bin/matlab"

clean-vectors source:
    rm -rf {{testdata}}/{{source}}

vectors-scipy: (clean-vectors "scipy")
    cd reference/scipy && uv run main.py {{testdata}}/scipy

vectors-matlab: (clean-vectors "matlab")
    {{matlab}} -sd reference/matlab -batch "generate('{{testdata}}/scipy', '{{testdata}}/matlab')"

vectors: vectors-scipy vectors-matlab

fixtures: vectors

lint:
    dprint check
    test -z "$(gofmt -l .)"
    go vet ./...
    cd {{scipy_reference_dir}} && uv run ruff format --check .
    cd {{scipy_reference_dir}} && uv run ruff check .

fmt:
    dprint fmt
    go fmt ./...
    cd {{scipy_reference_dir}} && uv run ruff format .
    cd {{scipy_reference_dir}} && uv run ruff check --fix .

test:
    go test ./...

test-verbose:
    go vet ./...
    go test ./... -v

tolerances:
    go run ./cmd/tolerances

bench:
    go test ./fft -run '^$' -bench . -benchmem -count 6 | tee temp/bench.txt

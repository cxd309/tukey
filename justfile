testdata := justfile_directory() / "testdata"
scipy_reference_dir := "reference/scipy"

clean-vectors source:
    rm -rf {{testdata}}/{{source}}

vectors-scipy: (clean-vectors "scipy")
    cd reference/scipy && uv run main.py {{testdata}}/scipy


vectors: vectors-scipy 

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

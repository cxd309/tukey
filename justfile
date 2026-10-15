testdata := justfile_directory() / "testdata"

clean-vectors source:
    rm -rf {{testdata}}/{{source}}

vectors-scipy: (clean-vectors "scipy")
    cd reference/scipy && uv run python -m vectors {{testdata}}/scipy

vectors: vectors-scipy 

fixtures: vectors

fmt:
    dprint fmt
    go fmt ./...
    cd testdata && uv run ruff format .

fixtures:
    cd testdata && uv run generate.py

test:
    go vet ./...
    go test ./...

test-verbose:
    go vet ./...
    go test ./... -v

tolerances:
    go run ./cmd/tolerances

bench:
    go test ./fft -run '^$' -bench . -benchmem -count 6 | tee temp/bench.txt
testdata := justfile_directory() / "testdata"
scipy_reference_dir := "reference/scipy"
# MATLAB on the PATH, or set MATLAB to its full path
matlab := env("MATLAB", "matlab")
# architectures to run on. each is skipped where this machine not able to run
goarchs := "amd64 arm64"

clean-vectors source:
    rm -rf {{testdata}}/{{source}}

vectors-scipy: (clean-vectors "scipy")
    cd reference/scipy && uv run main.py {{testdata}}/scipy

vectors-matlab: (clean-vectors "matlab")
    {{matlab}} -sd reference/matlab -batch "generate('{{testdata}}/scipy', '{{testdata}}/matlab')"

vectors: vectors-scipy vectors-matlab

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
# run command once per architecture in all goarchs that can be run on this machine
# result can differ between them (FMA, math assembly), so each is worth checking
# go vet first, so build errors are not mistaken for "can't run"
[private]
each-arch +command:
    #!/usr/bin/env bash
    set -euo pipefail
    for arch in {{goarchs}}; do
        export GOARCH=$arch
    go vet ./...
        if ! go test -run '^$' ./... >/dev/null 2>&1; then
            echo "== $arch: skipped, not able to be run on this machine"
            continue
        fi
        echo "== $arch"
        {{command}}
    done

test: (each-arch "go test ./...")

test-verbose: (each-arch "go test ./... -v")

tolerances: (each-arch "go run ./cmd/tolerances")

bench:
    go test ./fft -run '^$' -bench . -benchmem -count 6 | tee temp/bench.txt

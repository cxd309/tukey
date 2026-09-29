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

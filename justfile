fmt:
    dprint fmt
    go fmt ./...
    cd testdata && uv run ruff format .

fixtures:
    cd testdata && uv run generate.py

test: fixtures
    go vet ./...
    go test ./...

test-verbose: fixtures
    go vet ./...
    go test ./... -v

tolerances:
    go run ./cmd/tolerances

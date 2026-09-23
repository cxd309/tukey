fmt:
    dprint fmt
    go fmt ./...
    cd testdata && uv run ruff format .

fixtures:
    cd testdata && uv run generate.py

test: fixtures
    go test ./... -v
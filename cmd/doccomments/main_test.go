package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// allowed is a small snippet set for the tests: π and ².
var allowed = map[rune]bool{'π': true, '²': true}

func TestLineIssue(t *testing.T) {
	cases := []struct {
		name, line, want string
	}{
		{"plain prose", "// Butter designs a filter.", ""},
		{"a snippet character", "// The phase is π·x²", `character '·' isn't`},
		{"allowed characters", "// e^(-π) is π²", ""},
		{"an em dash", "// a filter — a design", `character '—' isn't`},
		{"a semicolon in prose", "// it rounds; it's exact", "semicolon in comment"},
		{"a semicolon in code", "//\tfor i := 0; i < n; i++ {", ""},
		{"an unpaired quote in code", "//\tx := 1 // the sample's value", "unpaired '"},
		{"paired quotes in code", "//\tifft(x, n, 'symmetric')", ""},
		{"an apostrophe in prose", "// the filter's gain", ""},
		{"an unpaired double quote in code", "//\ts := \"open", `unpaired "`},
	}
	for _, c := range cases {
		got := lineIssue(c.line, allowed)
		if c.want == "" {
			assert.Empty(t, got, c.name)
			continue
		}
		assert.True(t, strings.HasPrefix(got, c.want), "%s: got %q", c.name, got)
	}
}

// snippetCharacters reads every body from VS Code's JSON, trailing commas and all.
func TestSnippetCharacters(t *testing.T) {
	path := filepath.Join(t.TempDir(), "maths.code-snippets")
	require.NoError(t, os.WriteFile(path, []byte(`{
  "pi":      { "prefix": ";pi",      "body": "π" },
  "squared": { "prefix": ";squared", "body": "²" },
}`), 0o600))
	got, err := snippetCharacters(path)
	require.NoError(t, err)
	assert.Equal(t, allowed, got)

	require.NoError(t, os.WriteFile(path, []byte(`{}`), 0o600))
	_, err = snippetCharacters(path)
	assert.Error(t, err, "a file with no bodies")
}

// fileIssues reports a doc comment split from its declaration, and each bad line, with its line number.
func TestFileIssues(t *testing.T) {
	path := filepath.Join(t.TempDir(), "example.go")
	require.NoError(t, os.WriteFile(path, []byte(`package example

// Lost is a doc comment cut off by a blank line.

// Kept is attached; but this line has a semicolon.
func Kept() {}

func inside() {
	// comments inside a function are attached, and checked line by line
	_ = 1
}
`), 0o600))
	issues, err := fileIssues(path, allowed)
	require.NoError(t, err)
	require.Len(t, issues, 2)
	assert.Contains(t, issues[0], "example.go:3: attached to no declaration: Lost is")
	assert.Contains(t, issues[1], "example.go:5: semicolon in comment")
}

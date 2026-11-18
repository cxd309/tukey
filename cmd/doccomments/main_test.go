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

func TestMarkdownLineIssue(t *testing.T) {
	cases := []struct {
		name, line string
		fenced     bool
		want       string
	}{
		{"plain prose", "Every function returns an error.", false, ""},
		{"a curly apostrophe", "the filter’s gain", false, `character '’' isn't`},
		{"an em dash in a code block", "x := 1 // a — b", true, `character '—' isn't`},
		{"an allowed character", "Nyquist is π rad/sample.", false, ""},
		{"a semicolon in prose", "It rounds; it's exact.", false, "semicolon in prose"},
		{"a semicolon in inline code", "Write `for i := 0; i < n; i++` instead.", false, ""},
		{"a semicolon in a code block", "for i := 0; i < n; i++ {", true, ""},
	}
	for _, c := range cases {
		got := markdownLineIssue(c.line, c.fenced, allowed)
		if c.want == "" {
			assert.Empty(t, got, c.name)
			continue
		}
		assert.True(t, strings.HasPrefix(got, c.want), "%s: got %q", c.name, got)
	}
}

// markdownIssues tracks fenced code blocks, and reports each bad line with its line number.
func TestMarkdownIssues(t *testing.T) {
	path := filepath.Join(t.TempDir(), "doc.md")
	require.NoError(t, os.WriteFile(path, []byte("# Title\n"+
		"\n"+
		"Prose; with a semicolon.\n"+
		"\n"+
		"```go\n"+
		"for i := 0; i < n; i++ {}\n"+
		"```\n"+
		"\n"+
		"Back in prose; again.\n"), 0o600))
	issues, err := markdownIssues(path, allowed)
	require.NoError(t, err)
	require.Len(t, issues, 2)
	assert.Contains(t, issues[0], "doc.md:3: semicolon in prose")
	assert.Contains(t, issues[1], "doc.md:9: semicolon in prose")
}

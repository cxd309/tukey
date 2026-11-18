// Command doccomments fails on comments and docs that break the house style, so godoc, editors and GitHub show them as written.
// In Go comments it fails on:
//   - a top-level comment attached to no declaration, e.g. a doc comment split from its declaration by a blank line,
//     of which godoc shows only the part next to the declaration
//   - a character that isn't ASCII or in .vscode/maths.code-snippets, which includes em dashes
//   - a semicolon in a prose line, as plain English has no place for one
//   - an unpaired ', " or backtick on an indented line, which godoc shows as code,
//     and which an editor's hover highlights as an unclosed literal
//
// In markdown it fails on a character that isn't ASCII or in the snippets, anywhere,
// and on a semicolon in prose, outside code blocks and inline code.
//
// It reads every .go and .md file outside hidden directories, temp and testdata.
//
//	go run ./cmd/doccomments
package main

import (
	"fmt"
	"go/ast"
	"go/parser"
	"go/token"
	"io/fs"
	"os"
	"path/filepath"
	"regexp"
	"strings"
)

// snippetsPath is the file listing the characters comments may use beyond ASCII.
const snippetsPath = ".vscode/maths.code-snippets"

func main() {
	allowed, err := snippetCharacters(snippetsPath)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	found := 0
	err = filepath.WalkDir(".", func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() && path != "." && (strings.HasPrefix(d.Name(), ".") || d.Name() == "temp" || d.Name() == "testdata") {
			return filepath.SkipDir
		}
		var issues []string
		switch {
		case d.IsDir():
			return nil
		case strings.HasSuffix(path, ".go"):
			issues, err = fileIssues(path, allowed)
		case strings.HasSuffix(path, ".md"):
			issues, err = markdownIssues(path, allowed)
		default:
			return nil
		}
		if err != nil {
			return err
		}
		for _, issue := range issues {
			fmt.Println(issue)
		}
		found += len(issues)
		return nil
	})
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if found > 0 {
		fmt.Fprintf(os.Stderr, "%d comment or doc issue(s): see the writing rules in docs/development.md\n", found)
		os.Exit(1)
	}
}

// snippetCharacters is the set of characters in the snippets file's bodies.
// The file is VS Code's JSON with trailing commas, so the bodies are found by pattern rather than parsed.
func snippetCharacters(path string) (allowed map[rune]bool, err error) {
	text, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	allowed = map[rune]bool{}
	for _, m := range regexp.MustCompile(`"body":\s*"([^"]*)"`).FindAllSubmatch(text, -1) {
		for _, r := range string(m[1]) {
			allowed[r] = true
		}
	}
	if len(allowed) == 0 {
		return nil, fmt.Errorf("%s: no snippet bodies found", path)
	}
	return
}

// fileIssues lists every comment issue in the file at path, as "path:line: issue: text".
func fileIssues(path string, allowed map[rune]bool) (issues []string, err error) {
	fset := token.NewFileSet()
	f, err := parser.ParseFile(fset, path, nil, parser.ParseComments)
	if err != nil {
		return nil, err
	}
	report := func(pos token.Pos, issue, text string) {
		issues = append(issues, fmt.Sprintf("%s:%d: %s: %s", path, fset.Position(pos).Line, issue, strings.TrimSpace(text)))
	}
	for _, c := range floatingComments(f) {
		report(c.Pos(), "attached to no declaration", strings.SplitN(c.Text(), "\n", 2)[0])
	}
	for _, group := range f.Comments {
		for _, c := range group.List {
			for i, line := range strings.Split(c.Text, "\n") {
				if issue := lineIssue(line, allowed); issue != "" {
					report(c.Pos()+token.Pos(lineOffset(c.Text, i)), issue, line)
				}
			}
		}
	}
	return
}

// lineOffset is the byte offset of line i within text, for reporting the line it's on.
func lineOffset(text string, i int) (offset int) {
	for range i {
		offset += strings.IndexByte(text[offset:], '\n') + 1
	}
	return
}

// lineIssue is what's wrong with one line of a comment, or "" if nothing is.
// An indented line, with a tab after the //, is code to godoc.
func lineIssue(line string, allowed map[rune]bool) (issue string) {
	text := strings.TrimPrefix(strings.TrimSpace(line), "//")
	text = strings.TrimPrefix(text, "/*")
	code := strings.HasPrefix(text, "\t")
	for _, r := range text {
		if r > 127 && !allowed[r] {
			return fmt.Sprintf("character %q isn't ASCII or in %s", r, snippetsPath)
		}
	}
	if !code && strings.Contains(text, ";") {
		return "semicolon in comment"
	}
	if code {
		for _, q := range []string{"'", `"`, "`"} {
			if strings.Count(text, q)%2 == 1 {
				return fmt.Sprintf("unpaired %s in a code line", q)
			}
		}
	}
	return
}

// floatingComments lists the file's top-level comments that are attached to no declaration.
// The package's doc comment, comments above the package clause,
// and comments inside a declaration are all attached.
func floatingComments(f *ast.File) (floating []*ast.CommentGroup) {
	attached := map[*ast.CommentGroup]bool{f.Doc: true}
	var spans [][2]token.Pos
	for _, d := range f.Decls {
		switch d := d.(type) {
		case *ast.FuncDecl:
			attached[d.Doc] = true
		case *ast.GenDecl:
			attached[d.Doc] = true
		}
		spans = append(spans, [2]token.Pos{d.Pos(), d.End()})
	}
	for _, c := range f.Comments {
		if attached[c] || c.Pos() < f.Name.Pos() || inside(c, spans) {
			continue
		}
		floating = append(floating, c)
	}
	return
}

// inside reports whether c lies within one of the spans.
func inside(c *ast.CommentGroup, spans [][2]token.Pos) (ok bool) {
	for _, s := range spans {
		if c.Pos() > s[0] && c.End() < s[1] {
			return true
		}
	}
	return
}

// markdownIssues lists every issue in the markdown file at path, as "path:line: issue: text".
// Characters are checked on every line, and semicolons only in prose,
// outside fenced code blocks and inline code.
func markdownIssues(path string, allowed map[rune]bool) (issues []string, err error) {
	text, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	fenced := false
	for i, line := range strings.Split(string(text), "\n") {
		if strings.HasPrefix(strings.TrimSpace(line), "```") {
			fenced = !fenced
			continue
		}
		if issue := markdownLineIssue(line, fenced, allowed); issue != "" {
			issues = append(issues, fmt.Sprintf("%s:%d: %s: %s", path, i+1, issue, strings.TrimSpace(line)))
		}
	}
	return
}

// inlineCode matches a markdown code span, which may hold a semicolon, as code does.
var inlineCode = regexp.MustCompile("`[^`]*`")

// markdownLineIssue is what's wrong with one line of markdown, or "" if nothing is.
// fenced says whether the line is inside a fenced code block.
func markdownLineIssue(line string, fenced bool, allowed map[rune]bool) (issue string) {
	for _, r := range line {
		if r > 127 && !allowed[r] {
			return fmt.Sprintf("character %q isn't ASCII or in %s", r, snippetsPath)
		}
	}
	if !fenced && strings.Contains(inlineCode.ReplaceAllString(line, ""), ";") {
		return "semicolon in prose"
	}
	return
}

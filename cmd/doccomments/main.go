// Command doccomments fails if any top-level comment is attached to no declaration.
// A doc comment split from its declaration by a blank line is one:
// godoc shows only the part next to the declaration, and the rest is lost.
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
	"strings"
)

func main() {
	found := 0
	err := filepath.WalkDir(".", func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() && path != "." && (strings.HasPrefix(d.Name(), ".") || d.Name() == "temp" || d.Name() == "testdata") {
			return filepath.SkipDir
		}
		if d.IsDir() || !strings.HasSuffix(path, ".go") {
			return nil
		}
		floating, err := floatingComments(path)
		if err != nil {
			return err
		}
		for _, f := range floating {
			fmt.Println(f)
		}
		found += len(floating)
		return nil
	})
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if found > 0 {
		fmt.Fprintf(os.Stderr, "%d comment(s) attached to no declaration: join each to its declaration, or remove it\n", found)
		os.Exit(1)
	}
}

// floatingComments lists the file's top-level comments that are attached to no declaration,
// as "path:line: first line of the comment".
// The package's doc comment, comments above the package clause,
// and comments inside a declaration are all attached.
func floatingComments(path string) (floating []string, err error) {
	fset := token.NewFileSet()
	f, err := parser.ParseFile(fset, path, nil, parser.ParseComments)
	if err != nil {
		return nil, err
	}
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
		first := strings.SplitN(strings.TrimSpace(c.Text()), "\n", 2)[0]
		floating = append(floating, fmt.Sprintf("%s:%d: %s", path, fset.Position(c.Pos()).Line, first))
	}
	return
}

// inside reports whether c lies within one of the spans
func inside(c *ast.CommentGroup, spans [][2]token.Pos) bool {
	for _, s := range spans {
		if c.Pos() > s[0] && c.End() < s[1] {
			return true
		}
	}
	return false
}

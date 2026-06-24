;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("book" "12pt" "a4paper" "oneside" "spanish")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("babel" "es-tabla" "es-nodecimaldot" "spanish" "es-noquoting") ("mathtools" "") ("amssymb" "") ("amsthm" "") ("physics" "") ("bm" "") ("caption" "") ("graphicx" "") ("booktabs" "") ("array" "") ("tabularx" "") ("ragged2e" "") ("enumitem" "") ("xcolor" "") ("listings" "") ("inconsolata" "") ("mdframed" "") ("thmtools" "") ("hyperref" "colorlinks=true" "linkcolor=blue" "citecolor=red" "urlcolor=teal") ("cleveref" "nameinlink" "noabbrev") ("biblatex" "backend=biber" "style=apa" "sorting=nyt") ("geometry" "")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "config/preamble"
    "chapters/marco-teorico"
    "book"
    "bk12"))
 :latex)


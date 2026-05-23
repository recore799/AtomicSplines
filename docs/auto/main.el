;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("book" "12pt" "a4paper" "oneside" "spanish")))
   (TeX-run-style-hooks
    "latex2e"
    "config/preamble"
    "chapters/framework"
    "book"
    "bk12"))
 :latex)


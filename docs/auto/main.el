;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("book" "12pt" "a4paper" "oneside" "spanish")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("geometry" "") ("fontenc" "T1") ("babel" "es-tabla" "es-nodecimaldot" "spanish" "es-noquoting") ("mathtools" "") ("amssymb" "") ("amsthm" "") ("physics" "") ("bm" "") ("caption" "") ("graphicx" "") ("booktabs" "") ("array" "") ("tabularx" "") ("ragged2e" "") ("enumitem" "") ("siunitx" "") ("multirow" "") ("longtable" "") ("todonotes" "colorinlistoftodos" "textsize=tiny") ("xcolor" "") ("listings" "") ("inconsolata" "") ("mdframed" "") ("thmtools" "") ("hyperref" "colorlinks=true" "linkcolor=blue" "citecolor=red" "urlcolor=teal") ("cleveref" "nameinlink" "noabbrev") ("biblatex" "backend=biber" "style=apa" "sorting=nyt")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "config/preamble"
    "chapters/introduccion"
    "chapters/marco-teorico"
    "chapters/metodo"
    "chapters/resultados"
    "chapters/discusion"
    "chapters/conclusiones"
    "appendices/apendice_unidades"
    "appendices/apendice_tensores"
    "appendices/apendice_splines"
    "appendices/apendice_poisson"
    "appendices/apendice_metodos_auxiliares"
    "book"
    "bk12")
   (LaTeX-add-environments
    '("proof" LaTeX-env-args ["argument"] 0)))
 :latex)


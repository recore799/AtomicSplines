;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "preamble"
 (lambda ()
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
    "fontenc"
    "babel"
    "mathtools"
    "amssymb"
    "amsthm"
    "physics"
    "bm"
    "caption"
    "graphicx"
    "booktabs"
    "array"
    "tabularx"
    "ragged2e"
    "enumitem"
    "siunitx"
    "multirow"
    "longtable"
    "todonotes"
    "xcolor"
    "listings"
    "inconsolata"
    "mdframed"
    "thmtools"
    "hyperref"
    "cleveref"
    "biblatex")
   (TeX-add-symbols
    '("tabularxcolumn" 1))
   (LaTeX-add-environments
    '("proof" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-bibliographies
    "bib/bibliografia")
   (LaTeX-add-xcolor-definecolors
    "gruvbg"
    "gruvfg"
    "gruvred"
    "gruvgreen"
    "gruvgray"
    "gruvblue"
    "gruvyellow")
   (LaTeX-add-thmtools-declaretheoremstyles
    "miestilo-plano"
    "miestilo-caja")
   (LaTeX-add-thmtools-declaretheorems
    "theorem"
    "lemma"
    "proposition"
    "corollary"
    "definition"
    "example"))
 :latex)


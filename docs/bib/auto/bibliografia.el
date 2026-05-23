;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "bibliografia"
 (lambda ()
   (LaTeX-add-bibitems
    "fischer1977"
    "Fischer1990"
    "TingYun2001"
    "Garza2012"
    "Pauli1940"))
 '(or :bibtex :latex))


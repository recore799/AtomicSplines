;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "metodo"
 (lambda ()
   (LaTeX-add-labels
    "ch:metodo"
    "eq:poisson-radial"
    "eq:robin")
   (LaTeX-add-environments
    '("proof" LaTeX-env-args ["argument"] 0)))
 :latex)


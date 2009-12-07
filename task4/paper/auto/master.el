(TeX-add-style-hook "master"
 (lambda ()
    (TeX-add-symbols
     "ct")
    (TeX-run-style-hooks
     "tikz"
     "dot2texi"
     "autosize"
     "datetime"
     "latex2e"
     "memoir10"
     "memoir"
     "10pt"
     "oneside"
     "a4paper"
     "final"
     "english"
     "env/packages"
     "env/forloop"
     "env/languages"
     "env/graphics"
     "env/math"
     "env/preamble")))


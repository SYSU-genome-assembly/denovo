(TeX-add-style-hook "pop"
 (lambda ()
    (LaTeX-add-bibliographies
     "elex")
    (LaTeX-add-labels
     "tab:marks")
    (TeX-run-style-hooks
     "biblatex"
     "style=phys"
     "backend=biber"
     "xeCJK"
     ""
     "latex2e"
     "art10"
     "article"
     "a4paper")))


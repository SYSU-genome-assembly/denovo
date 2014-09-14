(TeX-add-style-hook "cum"
 (lambda ()
    (LaTeX-add-bibliographies
     "elex")
    (LaTeX-add-labels
     "t:notas"
     "fig:my_label")
    (TeX-run-style-hooks
     "geometry"
     "right=1in"
     "left=1in"
     "bottom=1.25in"
     "top=1.25in"
     "subcaption"
     "wrapfig"
     "graphicx"
     "balance"
     "amsmath"
     "xeCJK"
     ""
     "latex2e"
     "art12"
     "article"
     "a4paper"
     "12pt")))


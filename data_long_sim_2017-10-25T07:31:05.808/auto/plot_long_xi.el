(TeX-add-style-hook
 "plot_long_xi"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "pgfplots"
    "amsmath"
    "amsfonts"
    "amssymb")))


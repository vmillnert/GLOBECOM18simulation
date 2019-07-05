(TeX-add-style-hook
 "plot_long_Rpath_Umean"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "pgfplots"
    "pgfplotstable")
   (LaTeX-add-labels
    "p1"
    "p2"
    "p3")))


(TeX-add-style-hook
 "plot_lambda_e_a"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "pgfplots")
   (LaTeX-add-labels
    "e"
    "a"
    "u")))


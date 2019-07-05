(TeX-add-style-hook
 "plot_lambda_d_o"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "pgfplots")
   (LaTeX-add-labels
    "d"
    "a")))


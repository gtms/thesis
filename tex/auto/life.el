(TeX-add-style-hook "life"
 (lambda ()
    (LaTeX-add-environments
     "docspec")
    (TeX-run-style-hooks
     "booktabs"
     "fontspec"
     "latex2e"
     "tufte-book10"
     "tufte-book"
     "")))


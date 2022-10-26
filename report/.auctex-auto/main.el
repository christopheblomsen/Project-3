(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "10pt" "english")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("babel" "english") ("xcolor" "usenames" "dvipsnames" "svgnames" "table") ("hyperref" "colorlinks") ("standalone" "subpreambles=true") ("biblatex" "backend=biber" "style=chicago-authordate")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "inputenc"
    "babel"
    "amsmath"
    "graphicx"
    "varioref"
    "verbatim"
    "amsfonts"
    "geometry"
    "enumerate"
    "commath"
    "textcomp"
    "listings"
    "siunitx"
    "float"
    "xcolor"
    "hyperref"
    "import"
    "xifthen"
    "pdfpages"
    "transparent"
    "tikz"
    "pgfplots"
    "cancel"
    "standalone"
    "biblatex"
    "tocloft"
    "subfiles")
   (TeX-add-symbols
    '("uvec" 1)
    '("dd" 1)
    "doubleunderline")
   (LaTeX-add-bibliographies
    "../../../refs/refs")
   (LaTeX-add-listings-lstdefinestyles
    "mystyle")
   (LaTeX-add-xcolor-definecolors
    "mygreen"
    "mylilas"
    "codegreen"
    "codegray"
    "codepurple"
    "backcolour"
    "listinggray"
    "lbcolor"))
 :latex)


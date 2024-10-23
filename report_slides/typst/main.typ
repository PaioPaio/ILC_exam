#import "@preview/polylux:0.3.1": * // not necessary, but package for slides (already inside the template)
#import "@preview/grape-suite:1.0.0": slides //import template

#import slides: *

#set text(lang: "en")

#set math.vec(delim: "[")

// Template settings
#show: slides.with(
  no: none,
  series: [ILC Exam Report],
  title: [Implementation of f-ILC],

  author: "Lorenzo Paiola",
  show-semester: false,
  show-date: false,
  email: link("mailto:lorenzo.paiola@iit.it"),
  box-task-title: none,
)

#slide[
  = Overview

  - Implementation of "Iterative Learning in Functional Space for Non-Square Linear Systems" by _C. Della Santina_ and _F. Angelini_ @dellasantinaIterativeLearningFunctional2021.

  - Julia @Julia-2017 Code found at https://github.com/PaioPaio/ILC_exam
]

#slide[
  = What is Iterative Learning Control ?
]



#slide[
  = What is functional ILC ?
  #side-by-side[

  ][
    #image("images/blockschemes.png")]
]

#slide[
  = What is ILC ?
]

#slide[
  #set text(size: 0.7em)
  #bibliography("ManoloILC.bib", style: "nature")
]
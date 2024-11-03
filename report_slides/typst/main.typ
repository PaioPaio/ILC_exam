#import "@preview/polylux:0.3.1": * // not necessary, but package for slides (already inside the template)
#import "@preview/grape-suite:1.0.0": slides //import template

#import slides: *

#set text(lang: "en", font: "New Computer Modern")

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
  box-task-title: none
)

#show figure.caption: set text(size: 14pt)

#slide[
  = Overview
  \
  - Brief summary of Iterative Learning Control

  - Implementation of "Iterative Learning in Functional Space for Non-Square Linear Systems" by _C. Della Santina_ and _F. Angelini_ @dellasantinaIterativeLearningFunctional2021.

  - Julia @Julia-2017 Code found at https://github.com/PaioPaio/ILC_exam
]

#slide[
  = Iterative Learning Control in a nutshell

  Iterative Learning Control @bristowSurveyIterativeLearning2006 (ILC) generally concerns the control of a repeated task. It does so by:

  - Closing the loop in the *Iteration Domain* rather than directly time
  - Learning just the *Feed-Forward Input*

  #task()[
    === Remark
    ILC assumes that only the initial state is the same at each iteration, no assumptions are made about the terminal state.
  ]

]

#slide[

  ==== What is Iterative Learning Control ?
  
  In a nutshell ILC : 
  - 
  -
  - 

]

#slide[
  == System Set up
  #task[
    === LTI Continuous Time System

    #side-by-side()[
      $dot(x)_j = A x_j + B u_j, quad y_j=C x_j$
    ][
      with $x_j in RR^(n), u_j in RR^(l), y_j in RR^m$
    ]
  ]
  This system is:
  - *Iterated* and $j$ indicates the repetition index
  - Usually *non-square*, i.e. $l != m$, more interesting is the case where the system is underactuated $l<m$
  - *Sampled* only at a finite number of time instants ${T^1, dots , T^o}$
]

#slide[
  = Functional ILC
  #side-by-side[
    #set text(size: 0.9em)
    == fILC Structure
    - $vec(alpha^1_j, dots.v , alpha^(m o)_j)$ vector of weights updated at each iteration $j$
    - $l$ basis functions for each weight
    - Reference given at discrete set of sampled times ${T^1, dots, T^o}$, ($T^0=0$)
    - $L in RR^(m o times m o)$ learning matrix s.t. $rho(I - L H)<1$
  ][
    #set text(size: 1em)
    #image("images/blockschemes.png", fit: "cover")
  ]
]

#slide()[
  = Simulations
  #side-by-side()[  
    == Carts

    #image("images/MSD.png")

    Mass-spring-damper system actuated just on the last cart

  ][
    == Basketball in the wind

  #side-by-side()[
    #image("images/freethrow_wind.svg")
  ][
    Basket hoops in 
  ]
  ]

]

#slide()[
  ===  5 Carts
  #side-by-side()[
    #figure(
      image("images/MSD_5carts_30iters.svg", fit: "contain"),
      caption: "30 Iterations",
    )
  ][
    #figure(
      image("images/MSD_5carts_300iters.svg", fit: "contain"),
      caption: "300 Iterations - Not much changes",
    )
  ]
]

#slide()[
  === 8 Carts
  #side-by-side()[
    #figure(
      image("images/MSD_8carts_30iters.svg", fit: "contain"),
      caption: "30 Iterations - Can't perform to specification",
    )
  ][
    #figure(
      image("images/MSD_8carts_300iters.svg", fit: "contain"),
      caption: "300 Iterations - Better, still not perfect ",
    )
  ]
]

#slide()[
  ===  Basketball - *$d tilde(T) = 0.2s$*
  #side-by-side()[
    #figure(
      image("images/basket_8iters_0.2dT.svg", fit: "contain"),
      caption: "8 Iterations",
    )
  ][
    #figure(
      image("images/basket_100iters_0.2dT.svg", fit: "contain"),
      caption: "100 Iterations",
    )
  ]
]

#slide()[
  ===  Basketball - *$d tilde(T)$* = 0.5s
  #side-by-side()[
    #figure(
      image("images/basket_8iters_0.5dT.svg", fit: "contain"),
      caption: "8 Iterations",
    )
  ][
    #figure(
      image("images/basket_100iters_0.5dT.svg", fit: "contain"),
      caption: "100 Iterations",
    )
  ]
]


#slide[
  #set text(size: 0.7em)
  #bibliography("ManoloILC.bib", style: "nature")
]
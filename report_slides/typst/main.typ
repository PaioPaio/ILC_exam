#import "@preview/polylux:0.3.1": * // not necessary, but package for slides (already inside the template)
#import "@preview/grape-suite:1.0.0": slides //import template

#import slides: *

#set text(lang: "en")

#set math.vec(delim: "[")
#set math.mat(delim: "[")
#set math.equation(numbering: "(1)")

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

#show figure.caption: set text(size: 14pt)

#slide[
  = Overview
  \
  - Brief summary of Iterative Learning Control

  - Implementation of "Iterative Learning in Functional Space for Non-Square Linear Systems" by _C. Della Santina_ and _F. Angelini_ @dellasantinaIterativeLearningFunctional2021.

  - Julia @Julia-2017 Code found at https://github.com/PaioPaio/ILC_exam
]

#slide[
  = Iterative Learning Control in a Nutshell

  Iterative Learning Control @bristowSurveyIterativeLearning2006 (ILC) generally concerns the control of a repeated task. It does so by:

  - Closing the loop in the *Iteration Domain* rather than directly time
  - Learning just the *Feed-Forward Input*

  #task()[
    === Remark
    ILC assumes that only the initial state is the same at each iteration, no assumptions are made about the terminal state.
  ]

]

#slide[

  === Iterative Learning Control in a Nutshell
  == System Set up
  \

  #task[
    === LTI Continuous Time System

    #side-by-side()[
      $dot(x)_j = A x_j + B u_j, quad y_j=C x_j$
    ][
      with $x_j in RR^(n), u_j in RR^(l), y_j in RR^m$
    ]
  ]
  This system is *Iterated* and $j$ indicates the repetition index.
]

#slide()[
  === Iterative Learning Control in a Nutshell
  === Usual Policies

  #task[
    #align(center)[
      #grid(
        columns: 3, rows: 2, gutter: 10mm,
        [
          ==== P-type
          $u_(j+1)= u_j + L e_j$
        ],
        [
          ==== D-type
          $u_(j+1)= u_j + L dot(e)_j$
        ],
        [
          ==== I-type
          $u_(j+1)= u_j + L (e_(j+1) - e_j)$
        ],
        grid.cell(
          colspan: 2,
          align: center,
          [
            ==== General Rule
            $u_(j+1)= u_j + sum_(k=0)^(r) L_k e^((k))_(j)$
          ],
        ),
        grid.cell(
          colspan: 1,
          align: center,
          [
            \
            $e_j = overline(y) - y_j, quad$ $L_k in RR^(l times m)$
          ],
        )
      )
    ]
  ]
  Where $r$ is the relative degree of the system, $e^((k))_(j)$ is the $k$-th derivative. Notice that we are working with functions.

]


#slide[

  === Iterative Learning Control in a Nutshell

  == What is Hard ?

  - Usually the system is *non-square*, i.e. $l != m$, more interesting is the case where the system is underactuated $l<m$.
  - We *sample* $overline(y)$ only at a finite number of time instants ${T^1, dots , T^o}$ and so we have no access to its derivatives or even the true function $overline(y)$.
  - As specified before, we're working with functions and learning functions directly is hard. #text(15pt)[Moreover, we abuse notation and omit the proper definition of functions and just specify the output domain, eg $y(t): RR_+arrow RR^m$ becomes $y in RR^m$.]
]

#slide[
  = Functional ILC
  #task[
    The idea is to not learn $u_j$ directly at each iteration, but to learn some vector of weights $alpha_j$ that is multiplied to a library of base functions $pi_j$.
  ]
  What we want is
  $
    lim_(j arrow infinity) y_j (T^k) = overline(y)^k, quad forall k in 1 dots o
  $<objective>
  where $overline(y)^k$ is the sampled output at time $T^k$.
]

#slide()[
  === Functional ILC - Main Theorem
  #task[
    #set text(size: 0.8em)
    For a functional basis $pi = [pi^1 dots pi^o] in RR^(l times m o)$ and a control function constructed as
    $
      u_j (t) = pi alpha_j, quad
      alpha_(j+1) = alpha_j + L vec(overline(y)^1 - y_j (T^1), dots.v , overline(y)^o - y_j (T^o)) in RR^(m o)
    $
    @objective is achieved if $rho (I - L H) <1$, with
    $
      H = vec( integral_0^(T^1) C e^(A(T^1 - tau)) B pi (tau) d tau, dots.v, integral_0^(T^o) C e^(A(T^o - tau)) B pi (tau) d tau ) in RR^(m o times m o).
    $<Hmat>
  ]
]

#slide()[
  === Functional ILC - $pi$ choice
  #set text(size: 0.8em)
  Any choice of $pi$ for which $H$ is full rank is fine.

  The paper @dellasantinaIterativeLearningFunctional2021 gives a possible construction that satisfies this condition:
  $
    pi^i (t) = cases(B^top e^(A^top (T^i - t))C^top &"if" T^(i-1) <= t <= T^i",", 0 in RR^(l times m) &"otherwise.")
  $<pinormal>
  In the case in which the control can be applied only in a limited timespan, the basis
  $
    pi^i (
      t
    ) = cases(B^top e^(A^top (T^i - t))C^top &"if" T^(i-1) <= t <= T^(i-1) + d tilde(T), 0 in RR^(l times m) &"otherwise")
  $<piconstr>
  is appropriate.

  The learning matrix is instead set to $L = (H^top H + S)^(-1) H^top in RR^(m o times m o), quad S = I dot 10^(-2).$
]

#slide[
  === Functional ILC - Summary
  #side-by-side[
    #image("images/blockschemes.png", fit: "cover")
  ][
    #set text(size: 0.9em)
    == fILC Structure
    - $vec(alpha^1_j, dots.v , alpha^(m o)_j)$ vector of weights updated at each iteration $j$
    - $l$ basis functions for each weight
    - Reference given at discrete set of sampled times ${T^1, dots, T^o}$, ($T^0=0$)
    - $L in RR^(m o times m o)$ learning matrix s.t. $rho(I - L H)<1$
  ]
]

#slide()[
  = Examples
  #side-by-side()[
    === Carts

    #image("images/MSD.png")

    *Mass-Spring-Damper* system actuated by a force on the last cart. Uses @pinormal.

  ][
    === Basketball in the wind

    #side-by-side()[
      \
      Juggle a *Basketball* in place a few times and then try a free throw. Uses @piconstr.

    ][
      #image("images/freethrow_wind.svg")
    ]
    #text(10pt)[The details of the systems can be found directly in the paper @dellasantinaIterativeLearningFunctional2021.]
  ]
]

#slide()[
  === Results - MSD, 5 Carts
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
  === Results - MSD, 8 Carts
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
  === Results - Basketball, $d tilde(T)$ = 0.2s
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
  === Results - Basketball, $d tilde(T)$ = 0.5s
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

#slide()[
  = Remarks
  \
  - If the number of carts is $>8$ fILC struggles to converge, for higher numbers it outrights fails even after $tilde 10^3$ Iterations
  - The success of the scheme in the basketball example is dependent on the integrator used to simulate the system.
    - This is most likely due to the step-like nature of the constrained basis $pi_i$ rendering the linear ODE stiff
    - Strangely even integrators made for stiff ODEs fails if the method order is set too high or the timestep adaptive
    - Integrators that worked were: Heun, ROCK2, Ralston
]


#slide[
  #set text(size: 0.7em)
  #bibliography("ManoloILC.bib", style: "nature")
]
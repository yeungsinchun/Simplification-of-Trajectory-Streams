#set text(
  font: "Times New Roman",
  size: 12pt
)

#show heading.where(level: 1): set text(size: 1em)
#show heading.where(level: 2): set text(size: 0.8em)

#show heading: set block(above: 2em, below: 2em)

#set par(
  leading: 1.5em, 
  spacing: 1.5em,
  justify: true,
  first-line-indent: 2em,
)

#page(numbering: none)[
  #align(center + horizon)[
    #text(size: 24pt, weight: "bold")[Simplification of Trajectory Streams \ Progress Report]

    #v(2em)
    
    Yeung Sin Chun \
    scyeung\@connect.ust.hk \
    #datetime.today().display("[month repr:long] [day], [year]")
    
    #v(4em)
    
    #block(width: 80%)[
      #align(left)[
        #text(weight: "bold")[Abstract]
        
        #par(first-line-indent: 0em)[
          The ubiquitous use of GPS sensors has enabled real-time tracking of vehicles, which in turn enables the collection of massive trajectory data. Yet, for a massive stream, sending all vertices may be highly wasteful since only a small percentage of points on the trajectory are significant to maintain the shape of the original stream. While there are previous algorithms that does trajectory simplifications, few of them offers error guarantee using Fréchet distance. The aim of this project is to implement and benchmark a new streaming algorithm with a theorectical guarantee on Fréchet distance. We benchmark the algorithm in terms of the number of points in the simplified curve and the Fréchet distance achieved.
        ]
      ]
    ]
  ]
]

= Introduction
Trajectory stream simplification is a critical task in the era of ubiquitous GPS tracking. As vehicles, mobile devices, and sensors generate massive volumes of location data in real-time, transmitting every single data point becomes bandwidth-inefficient. Much of this data is redundant. For instance, a vehicle moving in a straight line generates many points that contribute little to the trajectory's overall shape. Simplification reduces this data volume while preserving the essential geometric features, enabling faster transmission, and more efficient real-time analytics.

While many software systems perform on-the-fly simplification, few algorithms offer rigorous quality guarantees that satisfy streaming requirements. In this project, we explore a streaming algorithm described in@algo. However, we will only study the algorithm in the context of $RR^2$ because of the complexity of implementing the algorithm in higher dimension.

For user-defined parameters $epsilon in (0, 1)$ and error bound $delta > 0$, the algorithm constructs a simplified curve $sigma$ in $RR^2$ that satisfies two key guarantees. First, the simplified curve is "close" to the original curve such that at any prefix of the original curve $tau[v_1, v_i]$, the simplified curve $sigma$ satisfies $d_F (sigma, tau[v_1,v_i]) <= (1 + epsilon)delta$. Second, the size of the simplified curve satisfies $|sigma| <= 2 dot "opt" - 2$ at any point during the algorithm, where $"opt"$ is the minimum number of vertices required to achieve a Fréchet error of at most $delta$ for the current prefix of the trajectory. The algorithm uses working storage of $O(epsilon^(-4))$ and each vertex in the original curve is processed in $O(epsilon^(-4)log 1/epsilon)$ time in $RR^2$.

== Fréchet Distance

The Fréchet distance is a measure of similarity between two curves that takes into account the location and ordering of the points along the curves. Let $S$ be a metric space. A curve $A$ in $S$ is a continuous map from the unit interval into $S$, i.e., $A: [0, 1] -> S$. A reparameterization $alpha$ of $[0, 1]$ is a continuous, non-decreasing, surjection $alpha: [0, 1] -> [0, 1]$.

Let $A$ and $B$ be two continuous curves in $S$. The Fréchet distance $d_F(A, B)$ is defined as the infimum over all reparameterizations $alpha$ and $beta$ of $[0, 1]$ of the maximum distance between $A(alpha(t))$ and $B(beta(t))$ for $t in [0, 1]$. Formally:

$ d_F(A, B) = inf_(alpha, beta) max_(t in [0, 1]) d(A(alpha(t)), B(beta(t))) $

where $d$ is the distance metric in $S$. We adopt the usual Euclidean distance.

Intuitively, this metric is often illustrated using the "dog-walking" analogy: imagine a person walking along curve $A$ and a dog walking along curve $B$. Both can control their speed but cannot move backwards. The Fréchet distance corresponds to the minimum length of the leash required to connect the dog and the person throughout their entire walk.

We will use the implementation in @frechet for measuring Fréchet distance.

= Other works
T
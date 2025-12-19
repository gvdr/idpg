#import "@preview/unequivocal-ams:0.1.2": ams-article, theorem, proof
#import "@preview/marginalia:0.2.3" as marginalia: note, notefigure
#show: marginalia.setup.with(
  inner: ( far: 5mm, width: 35mm, sep: 1mm ),
  outer: ( far: 5mm, width: 35mm, sep: 1mm )
)

#let mred(x) = text(fill: red)[$#x$]
#let mgreen(x) = text(fill: green)[$#x$]
#let mblue(x) = text(fill: blue)[$#x$]
#let mbred(x) = math.bold([$#text(fill: red)[$#x$]$])
#let mbgreen(x) = math.bold([$#text(fill: green)[$#x$]$])
#let mbblue(x) = math.bold([$#text(fill: blue)[$#x$]$])



#let date = datetime.today().display("[month repr:long] [day], [year]")

#show: ams-article.with(
  title: [Intensity Dot Product Graphs],
  authors: (
    (
      name: "Giulio Valentino Dalla Riva",
      organization: [Baffelan OU],
      location: [Noumea, New Caledonia],
      email: "me@gvdallariva.net",
      url: "www.baffelan.com"
    ),
    (
      name: "Matteo dalla Riva",
      department: [Department of Mathematics],
      organization: [Università di Padova],
      location: [Vicenza]
    ),
  ),
  abstract: [Random graphs have been successfully used to model a large variety of natural, cultural, social, economic, and other phenomena. In this work we build on the Random Dot Product Graph Models (RDPG), which have been well studied in the statistical literature. In an RDPG model, while the presence of a link between each pair of nodes in a graph is governed by a certain probability distribution, nodes are given as essentially deterministic objects. This is the case in all other random graph models studied so far. In this manuscript we study a model for random graphs in which nodes, and their propensity to form links, are governed by an intensity measure defined on a continuous Euclidean space. We also extend this framework to temporal random graphs, and show how to describe the temporal evolution of graphs in this model in terms of classic partial differential equations regimes.],
)

= A summary of the idea

A graph is a collection of some items that have (or not) established a pairwise connection (@graph). We call the items nodes, the connections edges, and we represent it as a set of dots linked by arrows.

Random graph models consider graphs as the outcome of a probability process @newman2018networks. In particular, most (or all) existing models for random graphs assume that nodes are somehow given, while an edge is established following a certain probability rule.

Here we extend the randomness to be inherent to the nodes as well. In particular, rather than assuming that nodes are given a priori, we consider nodes as points sampled from a Poisson point process, with a certain intensity, somewhere in a certain continuous, metric space. Then, the probability of establishing a connection is determined by the position of the sample points in that space.

In so doing, both the number and identity of the points are not determined _a priori_, but rather the outcome of a stochastic process.

Moreover, we can use classic differential equation techniques to study the time (or space) evolution of the intensities, and see how that affects the observed graphs.

= Random Dot Product Graphs

Let $G=(V,E)$ be a simple, directed graph, with nodes $i,j in V = { 1, 2, … } $ and edges $(i arrow j) in E subset V times V$ (we admit self connecting edges, i.e., $i arrow i$ is a possible edge). #notefigure(
  image("Graph.png", width: 200%),
  label: <graph>,
  caption: [A simple, directed graph with 5 nodes and a bunch of edges.],
)

We'll consider graphs as the outcome of random processes.
This means that we associate to any possible graph a certain probability of being observed.
In particular, let $V$ be a given set ${1, 2, dots, N}$ of nodes, every ordered couple in $V times V$ is in $E$ with a certain probability $p_"ij"$; we define the matrix of interaction probabilities $bold(P)$ so that $bold(P)_"ij" = p_"ij"$, and denote it $bold(P)_G$ if we need to be explicit regarding what graph it is associated with.

Notice here that $bold(P)$ completely determines the probability of observing a given graph $G=(V,E)$: the probability of $G$ will be given by the probability of observing exactly the links in $E$ and not observing the links not in $E$.

Random graph models are described by how they determine those interaction probabilities.

== RDPG as generating model

In a Random Dot Product Graph (RDPG) model @young2007random, each node is associated with two $d$-dimensional vectors, $mgreen(arrow(g_i)) = g(i)$ and $mred(arrow(r_i)) = r(i)$ (@node_vectors). The vectors $g(dot)$ and $r(dot)$ are chosen so that $g(i) dot r(j) in [0,1]$ for every $(i,j) in V times V$. A pair of nodes $i,j$ is an edge in $E$ with probability $p_(i,j) = mgreen(arrow(g_i))  dot mred(arrow(r_j)) = g(i) dot r(j)$.
#notefigure(
  image("node_vectors.png", width: 200%),
  label: <node_vectors>,
  caption: [Each node of a graph is associated with a pair of vectors, one green and one red. The probability of observing an edge between two nodes is given by the dot product of a green and a red vector.],
)

We can then consider two matrices $mbgreen(G)$ and $mbred(R)$, where the rows $mbgreen(G)_(i,dot)$ of $mbgreen(G)$ are the vectors $g(i)$ and the columns $mbred(R)_(dot,i)$ of $mbred(R)$ are the vectors $r(i)$ for every $i$ in $V$. We have that the matrix multiplication
$
mbgreen(G) mbred(R) = bold(P)
$
and, hence, the two matrices $(mbgreen(G), mbred(R))$ contain all the information of the random graph model (the number of the nodes is given by the number of rows of $mbgreen(G)$, that is the number of columns of $mbred(R)$).


It is convenient for our intuition to consider the vectors $mgreen(arrow(g_i))$ and $mred(arrow(r_i))$ as the node $i$ propensity to either propose or accept an edge connection (@points).
Furthermore, it is convenient to visualize each node $i$ as a pair of points in two $d$ dimensional metric spaces, that we will refer to (with some abuse of notation) as the green space $mgreen(G)$ and the red space $mred(R)$.
The coordinate of $i$ in these two spaces is given by $mgreen(arrow(g_i))$ and $mred(arrow(r_i))$. Thus, we can read the function $g(dot)$ and $r(dot)$ as the projection of the nodes to $mgreen(G)$ and $mred(R)$ respectively. Hence, we can also see an RDPG model $(mbgreen(G), mbred(R))$ as defined by a given set of $N$ points in the spaces $mgreen(G)$ and $mred(R)$, which offers a nice geometric representation of the random graph model. #figure(
  image("points.png", width: 70%),
  caption: [The green and red vectors associated with the nodes define a set of points in two metric spaces. We define two matrices: one green, whose rows will be the coordinates of the green points, one red, whose columns will be the coordinates of the red points. Their matrix multiplication gives all the connection probabilities.],
)<points>

To summarize, under a RDPG model, nodes are associated with points in a certain pair of spaces $mgreen(G)$ and $mred(R)$, and the probability of observing an edge between two points is given by the dot product $mgreen(arrow(g_i)) dot mred(arrow(r_j))$.

== Inference of an RDPG model

The inference of RDPG model parameters goes the other way round than the generation task.

Given an observed graph $G=(V,E)$ we are posed with the problem of identifying the two most likely matrices $(mbgreen(G), mbred(R))$ that generated $G$.

This will be accomplished in two steps:
1. Infer the right dimension $d$ for the spaces $mgreen(G)$ and $mred(R)$ (notice: the dimension, not the number of points)
2. Infer the positions of the $N$ points in $mgreen(G)$ and $mred(R)$.

First of all, we need to define the adjacent matrix of $G$. This is matrix $bold(A)$ such that
$
bold(A)_(i,j) = cases(
  1 "if" i arrow j in E,
  0 "otherwise" #h(1cm) .
) 
$

Classic results @athreya2018statistical show that both (i) and (ii) can be achieved with elementary linear algebra tools, namely some small manipulation of the singular value decomposition of $bold(A)$ (and indeed (i) is the most challenging!)

Let $bold(A) = bold(U) bold(Sigma) bold(V)^T$ be the singular value decomposition of $bold(A)$, that is $bold(U)$ and $bold(V)$ are orthogonal matrices and $bold(Sigma)$ is an $N times N$ diagonal matrix with non-negative real coefficients on its diagonal in decreasing order#footnote[The singular value decomposition is not unique: any orthogonal transformation of $bold(U)$ and $bold(V)$ would provide an equivalent decomposition; the singular values $bold(Sigma)$, on the other hand, are uniquely determined.]. The elements $sigma_i = bold(Sigma)_(i,i)$ are known as the singular values of $bold(A)$.

An optimal dimension $hat(d)$ can be inferred solely from the sequence of singular values $sigma_i$. There are various techniques for doing it, and the technicality is left to the curious reader (see for example @zhu2006automatic @gavish2014optimal @chatterjee2015matrix).

Let's define the two matrices $mbgreen(tilde(G)) = bold(U)|_hat(d) sqrt(bold(Sigma)|_hat(d))$ and $mbred(tilde(R)) = sqrt(bold(Sigma)|_hat(d)) (bold(V)|_hat(d))^T$, where $M|_k$ is the truncation of a matrix $M$ to its first $k$ columns, and $sqrt(bold(Sigma))_(i,i) = sqrt(sigma_i)$ is the element-wise square root of $bold(Sigma)$.

Then, we have that
$
cases(
  mbgreen(G) tilde.eq mbgreen(tilde(G)),
  mbred(R) tilde.eq mbred(tilde(R)) #h(1cm) .
)
$
In particular, the matrix $bold(tilde(A)) = mbgreen(tilde(G)) mbred(tilde(R)) tilde.eq bold(A)$
is optimal in the sense that it minimizes the distance to $bold(A)$ in Frobenius norm, that is:
$
bold(tilde(A)) = limits(arg min)_(bold(Mu) "of rank" hat(d)) norm(bold(A) - bold(Mu))_F #h(0.3cm) .
$

To summarize, given an observed graph $G=(V,E)$, a little bit of linear algebra, namely a singular value decomposition, is all that is needed to infer the parameters $(mbgreen(G), mbred(R))$ of the most likely RDPG model to generate $G$.

= Intensities

Notice that in a RDPG model, while the edges are probabilistic, the nodes are not: their number and their identities, that is their propensities to propose and accept an edge, are completely determined by the model parameters#footnote[And the same can be said for every other random network model.].

Now, we move a step forward in the probabilistic modelling of a network, dropping that deterministic nature of nodes.

To do this, we introduce the novel family of #emph[Intensity Dot Product Graph] (IDPG) models.

== The latent space

Before defining an IDPG, we address a technical constraint. In an RDPG, the vectors $g(i)$ and $r(j)$ must satisfy $g(i) dot r(j) in [0,1]$ for all pairs of nodes, so that the dot product can be interpreted as an edge probability. 

We ensure this by restricting the latent spaces $mgreen(G)$ and $mred(R)$ to be subsets of the $(d-1)$-dimensional probability simplex:
$ Delta^(d-1) = { x in bb(R)^d : x_k >= 0 "for all" k, sum_(k=1)^d x_k = 1 } $

For any $mgreen(arrow(g)), mred(arrow(r)) in Delta^(d-1)$, we have $mgreen(arrow(g)) dot mred(arrow(r)) in [0,1]$, as required.

Since all observable quantities in the model depend only on inner products between vectors, the model is invariant under orthogonal transformations of the latent space. The simplex provides a canonical representation that also guarantees valid connection probabilities.

#note[The coordinates on the simplex admit a natural interpretation as mixed membership weights: a node $i$ belongs to latent "community" $k$ with weight $g_(i,k)$. This connects IDPG to mixed-membership stochastic blockmodels.]

== IDPG as generating model

In an RDPG model, nodes are given as a fixed set of points in $mgreen(G)$ and $mred(R)$. In an IDPG model, we replace this discrete set with a continuous intensity measure, and nodes (or interaction sites) emerge as samples from a Poisson point process.

=== The intensity measure

Let $Omega = mgreen(G) times mred(R) subset.eq Delta^(d-1) times Delta^(d-1)$ be the product of the two latent spaces. An IDPG is specified by an #emph[intensity measure] $lambda$ on $Omega$. We assume $lambda$ admits a density $mbblue(rho)$ with respect to Lebesgue measure, so that for any region $A subset.eq Omega$:
$ lambda(A) = integral_A mbblue(rho)(mgreen(arrow(g)), mred(arrow(r))) dif mgreen(arrow(g)) dif mred(arrow(r)) $

Crucially, $mbblue(rho)$ is not a probability density: it need not integrate to one. The total intensity $lambda(Omega) = integral_Omega mbblue(rho)$ gives the expected number of sampled points.

=== Interaction sites and entity persistence

A sample from the Poisson point process PPP($lambda$) yields a random collection of points in $Omega$. Each sampled point $(mgreen(arrow(g)), mred(arrow(r)))$ represents an #emph[interaction site]: an entity that can propose connections via its green coordinate $mgreen(arrow(g))$ and accept connections via its red coordinate $mred(arrow(r))$.

The structure of the resulting graph depends on how long these entities #emph[persist] — that is, how many interactions each entity participates in before disappearing. This persistence parameter interpolates between two limiting regimes:

#block(
  fill: luma(245),
  inset: 10pt,
  radius: 4pt,
)[
*Definition (IDPG).* Let $lambda$ be a finite intensity measure on $Omega = mgreen(G) times mred(R) subset.eq Delta^(d-1) times Delta^(d-1)$, with density $mbblue(rho)$. An Intensity Dot Product Graph with intensity $lambda$ is generated by sampling interaction sites from PPP($lambda$). The probability that an entity at $(mgreen(arrow(g)), mred(arrow(r)))$ forms a directed edge to an entity at $(mgreen(arrow(g')), mred(arrow(r')))$ is given by the dot product $mgreen(arrow(g)) dot mred(arrow(r'))$.

The realized graph depends on entity persistence:
- *Long-lived entities*: sampled points persist, and all pairs of entities may interact.
- *Ephemeral entities*: each sampled point exists only for a single interaction.
]

These two regimes yield different graph structures, which we now characterize.

=== Long-lived entities (node-centric view)

In this regime, sampled points persist as #emph[nodes]. Each node $i$ is characterized by its coordinates $(mgreen(arrow(g)_i), mred(arrow(r)_i))$ in the latent space. Given $N$ sampled nodes, every ordered pair $(i, j)$ may form an edge $i arrow j$ independently with probability $mgreen(arrow(g)_i) dot mred(arrow(r)_j)$.

This recovers a classical graph structure: nodes can participate in multiple edges (as source via $mgreen(arrow(g))$, as target via $mred(arrow(r))$), giving rise to nontrivial degree distributions, paths, triangles, and connected components.

#figure(
  image("densities.png", width: 100%),
  caption: [From a set of points, we move toward intensity functions defining the expected density of interaction sites in the latent space. Regions of higher intensity will contain more sites on average. The green coordinate describes propensity to propose connections; the red coordinate describes propensity to accept them.],
)<density>

The expected number of nodes is:
$ bb(E)[N] = lambda(Omega) = integral_Omega mbblue(rho)(mgreen(arrow(g)), mred(arrow(r))) dif mgreen(arrow(g)) dif mred(arrow(r)) $

The expected number of edges involves a double integral over pairs of interaction sites:
$ bb(E)[|E|] = integral_Omega integral_Omega (mgreen(arrow(g)) dot mred(arrow(r'))) lambda(dif mgreen(arrow(g)), dif mred(arrow(r))) lambda(dif mgreen(arrow(g')), dif mred(arrow(r'))) $
where the connection probability $mgreen(arrow(g)) dot mred(arrow(r'))$ depends on the green coordinate of the source and the red coordinate of the target.

#note[The intensity $lambda$ introduces a natural scale: multiplying $lambda$ by a constant $c$ scales E[$N$] by $c$ and E[$|E|$] by $c^2$.]

=== Ephemeral entities (edge-centric view)

In this regime, each sampled point exists only long enough to form a single interaction. A sampled point $(mgreen(arrow(g)), mred(arrow(r')))$ is interpreted directly as an #emph[edge]: a connection from a source "at $mgreen(arrow(g))$" to a target "at $mred(arrow(r'))$".

This can be understood as sampling two ephemeral nodes — a source contributing only its green coordinate $mgreen(arrow(g))$, and a target contributing only its red coordinate $mred(arrow(r'))$. The source's red coordinate and the target's green coordinate are never used, as these entities disappear immediately after their single interaction.

Under this interpretation, the intensity $mbblue(rho)(mgreen(arrow(g)), mred(arrow(r)))$ represents the #emph[opportunity structure] for interactions: the rate at which potential interactions between sources at $mgreen(arrow(g))$ and targets at $mred(arrow(r))$ arise. This opportunity structure need not factorize — correlations in $mbblue(rho)$ can capture environmental, spatial, or social constraints that make certain (source, target) pairings more or less likely to encounter each other. For instance, at a social gathering, who is present, who occupies which part of the room, and who arrived together all create correlations in which pairs of individuals have the opportunity to interact.

The edge intensity — accounting for the connection probability once an opportunity arises — is:
$ mbblue(rho_E)(mgreen(arrow(g)), mred(arrow(r))) = mbblue(rho)(mgreen(arrow(g)), mred(arrow(r))) dot.c (mgreen(arrow(g)) dot mred(arrow(r))) $

This factorization reveals a two-stage sampling interpretation: first, potential interactions are drawn from the opportunity measure $mbblue(rho)$; second, each opportunity $(mgreen(arrow(g)), mred(arrow(r)))$ converts to an actual edge with probability $mgreen(arrow(g)) dot mred(arrow(r))$. Sampling directly from $mbblue(rho_E)$ is equivalent to this two-stage process, but the decomposition clarifies the distinct roles of opportunity structure and connection probability.

An edge-centric IDPG is then a Poisson point process on $mgreen(G) times mred(R)$ with intensity $mbblue(rho_E)$.

The expected number of edges is:
$ bb(E)[L] = integral_Omega mbblue(rho_E)(mgreen(arrow(g)), mred(arrow(r))) dif mgreen(arrow(g)) dif mred(arrow(r)) = integral_Omega mbblue(rho)(mgreen(arrow(g)), mred(arrow(r))) dot.c (mgreen(arrow(g)) dot mred(arrow(r))) dif mgreen(arrow(g)) dif mred(arrow(r)) $

#note[In the edge-centric view, the probability of two edges sharing an endpoint is zero (in continuous space). The resulting structure is a collection of disconnected edges rather than a connected graph. Traditional graph structure can be recovered through discretization or clustering of the latent space.]

=== Relating the two views

The node-centric and edge-centric views share the same mathematical specification — an intensity $mbblue(rho)$ on $mgreen(G) times mred(R)$ — but differ in interpretation:

#block(
  fill: luma(250),
  inset: 8pt,
  radius: 4pt,
)[
#table(
  columns: (auto, 1fr, 1fr),
  inset: 8pt,
  align: left,
  [], [*Node-centric*], [*Edge-centric*],
  [Sampled point $(mgreen(arrow(g)), mred(arrow(r)))$], [A node with coordinates $(mgreen(arrow(g)), mred(arrow(r)))$], [An edge from $mgreen(arrow(g))$ to $mred(arrow(r))$],
  [$mbblue(rho)(mgreen(arrow(g)), mred(arrow(r)))$], [Intensity of nodes at $(mgreen(arrow(g)), mred(arrow(r)))$], [Opportunity intensity for edges $mgreen(arrow(g)) arrow mred(arrow(r))$],
  [Correlations in $mbblue(rho)$], [A node's proposing and accepting roles are related], [Certain (source, target) pairs more likely to meet],
  [Entity persistence], [Long-lived (many interactions)], [Ephemeral (one interaction)],
  [Scaling of E[edges]], [Quadratic in intensity], [Linear in intensity],
)
]

These are limiting cases of entity persistence. Writing $tau$ for the lifetime of an entity (the duration over which it can participate in interactions):
- *Edge-centric*: entities persist for an instant ($tau arrow 0$)
- *Node-centric*: entities persist indefinitely ($tau arrow infinity$)

An intermediate regime — where entities form $k$ connections before disappearing — interpolates between these limits. This perspective connects naturally to temporal dynamics: if the intensity $mbblue(rho)$ evolves in time, the appropriate view depends on the timescale of entity persistence relative to the timescale of intensity evolution.

=== The product case

A case of particular interest arises when the intensity measure factorizes:
$ mbblue(rho)(mgreen(arrow(g)), mred(arrow(r))) = mbgreen(rho_G)(mgreen(arrow(g))) dot.c mbred(rho_R)(mred(arrow(r))) $

The product form has a natural interpretation in each view:
- *Node-centric*: a node's propensity to propose connections (its position in $mgreen(G)$) is independent of its propensity to accept connections (its position in $mred(R)$).
- *Edge-centric*: sources and targets encounter each other uniformly at random — there is no preferential opportunity structure pairing certain sources with certain targets.

Let us define:
- the marginal total intensities $mgreen(c_G) = integral_mgreen(G) mbgreen(rho_G)(mgreen(arrow(g))) dif mgreen(arrow(g))$ and $mred(c_R) = integral_mred(R) mbred(rho_R)(mred(arrow(r))) dif mred(arrow(r))$,
- the intensity-weighted mean positions $mbgreen(mu_G) = integral_mgreen(G) mgreen(arrow(g)) mbgreen(rho_G)(mgreen(arrow(g))) dif mgreen(arrow(g))$ and $mbred(mu_R) = integral_mred(R) mred(arrow(r)) mbred(rho_R)(mred(arrow(r))) dif mred(arrow(r))$.

==== Node-centric (long-lived)

The expected number of nodes is $bb(E)[N] = mgreen(c_G) dot.c mred(c_R)$.

Defining the normalized mean positions $mbgreen(tilde(mu)_G) = mbgreen(mu_G) \/ mgreen(c_G)$ and $mbred(tilde(mu)_R) = mbred(mu_R) \/ mred(c_R)$, the expected edge count is (see @appendix:derivations):
$ bb(E)[|E|] = (bb(E)[N])^2 dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) $
The expected number of edges is the square of the expected number of nodes, scaled by an average connection probability $mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R) in [0,1]$.

==== Edge-centric (ephemeral)

In the product case, the edge intensity becomes:
$ mbblue(rho_E)(mgreen(arrow(g)), mred(arrow(r))) = mbgreen(rho_G)(mgreen(arrow(g))) dot.c mbred(rho_R)(mred(arrow(r))) dot.c (mgreen(arrow(g)) dot mred(arrow(r))) $

The expected number of edges is (see @appendix:derivations):
$ bb(E)[L] = mbgreen(mu_G) dot mbred(mu_R) = mgreen(c_G) dot.c mred(c_R) dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) $

The quantity $mgreen(c_G) dot.c mred(c_R)$ appears in both views but with different interpretations:
- *Node-centric*: the expected number of persistent nodes, $bb(E)[N]$.
- *Edge-centric*: the expected number of source-target #emph[opportunities] — encounters between ephemeral sources and targets, before applying the connection probability.

Since both views share the same underlying intensity $mbblue(rho)$, we can use $bb(E)[N] = mgreen(c_G) dot.c mred(c_R)$ as a common scale parameter, giving:

$ bb(E)[|E|] &= (bb(E)[N])^2 dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) #h(1cm) &"(node-centric)" \
bb(E)[L] &= bb(E)[N] dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) #h(1cm) &"(edge-centric)" $

The ratio $bb(E)[|E|] \/ bb(E)[L] = bb(E)[N]$ captures the key distinction: in the node-centric view, each of the $bb(E)[N]$ nodes participates in $O(N)$ potential edges (as both source and target), yielding quadratic scaling. In the edge-centric view, each of the $bb(E)[N]$ opportunities yields at most one edge, giving linear scaling.

== Inference of an IDPG model

The inference problem for IDPG is: given an observed graph $G = (V, E)$, can we recover the intensity $mbblue(rho)$, or at least the marginal intensities $mbgreen(rho_G)$ and $mbred(rho_R)$ under a product assumption?

For a node-centric IDPG, a natural approach proceeds in two stages:

+ *Embed the observed nodes.* Apply the standard RDPG inference procedure: compute the singular value decomposition of the adjacency matrix $bold(A)$, select an appropriate dimension $hat(d)$, and obtain estimated positions $(hat(mgreen(arrow(g)_i)), hat(mred(arrow(r)_i)))$ for each observed node.

+ *Estimate the intensity.* Treat the embedded positions as a point cloud and apply density estimation techniques @silverman2018density to recover $mbblue(rho)$, or its marginals $mbgreen(rho_G)$ and $mbred(rho_R)$ under a product assumption.

The feasibility and accuracy of this procedure may depend on additional assumptions about the structure of $mbblue(rho)$. A natural constraint is to model $mbblue(rho)$ as a mixture of multivariate Gaussian distributions, which offers a flexible yet tractable family for density estimation.

For edge-centric observations, the inference problem requires specifying the observation model. If edges are observed directly as pairs of latent positions, estimating $mbblue(rho_E)$ reduces to density estimation on the point cloud. More commonly, one observes a discretized or aggregated version — for instance, interaction counts between clusters or categories — requiring additional modeling to relate the summary to the underlying continuous intensity.

= Time indexing and beyond

So far we considered graphs, and their random models, as stuck in time.

However, many applications motivate the study of graphs as dynamical processes. This means that instead of observing one graph $G$, we observe a sequence of graphs $G_t = (V_t,E_t)$, where the index $t$ is commonly assumed to stand for time. In most RDPG time extensions, we do consider $V_t = V$ for each time $t$, that is, the set of vertices does not change, while the connections between them can change.

Each graph $G_t$ is generated by an RDPG model with parameters  $mbgreen(G)_t$ and $mbred(R)_t$. A body of statistical results guides us to decide whether, given two graphs $G_t$ and $G_(t + delta t)$, we can determine that $(mbgreen(G),mbred(R))_t = (mbgreen(G),mbred(R))_(t + delta t)$ or not, that is, whether the change in the observation is actually induced by a movement of the points in $mgreen(G)$ and $mred(R)$ or by the inherent variability of the observation process.

Eventually, but this has so far received less attention, the graphs can be indexed by more than one variable, e.g., we could consider a spatiotemporal distribution of graphs $G_(x,y,t)$ where $x,y$ are some geographic coordinates and $t$ is time.

So far, there hasn't been an attempt to study from a dynamical system perspective the movement of the graph points in time and space.

== Two temporal scales

Introducing time into IDPG reveals two distinct dynamical scales:

=== Sampling dynamics (fast scale)

The Poisson point process describes when and where interaction sites appear. PPP is memoryless: events occur independently at rate $mbblue(rho)$. But richer temporal structure can arise from generalizing the sampling process:

- *Hawkes processes*: self-exciting point processes where past events increase the rate of future events. An interaction at $(mgreen(arrow(g)), mred(arrow(r)))$ temporarily boosts the intensity nearby, producing temporal clustering. This could model social reinforcement or predator-prey encounter dynamics.

- *Cox processes*: doubly stochastic processes where the intensity $mbblue(rho)$ is itself random, introducing additional variability in event rates.

These generalizations govern the temporal correlations in _when_ entities appear, while the spatial structure of $mbblue(rho)$ governs _where_ they appear.

=== Intensity evolution (slow scale)

On a slower timescale, the intensity $mbblue(rho)$ itself can evolve. We write $mbblue(rho)(mgreen(arrow(g)), mred(arrow(r)), t)$ for a time-varying intensity. The sampling process then operates on a landscape that drifts over time.

To handle the simplex constraint, we embed $Delta^(d-1)$ in $bb(R)^d$ and use standard PDE operators. The boundary of the simplex requires care. Two natural choices arise:

- *Absorbing boundary* ($mbblue(rho) = 0$ on $partial Delta$): intensity that reaches the boundary vanishes. Total mass decreases over time, and $bb(E)[N(t)]$ shrinks. This models extinction or loss of entities that become too extreme in their interaction propensities.

- *Reflecting boundary* (no-flux condition $nabla mbblue(rho) dot arrow(n) = 0$): intensity cannot escape the simplex. Total mass is conserved, so $bb(E)[N(t)]$ remains constant even as the distribution evolves. This may be more natural for ecological applications where species cannot leave the space of viable niches — they accumulate at boundaries rather than disappearing.

In either case, no interaction sites are ever sampled at positions yielding invalid connection probabilities (outside the simplex, the intensity is zero or inaccessible).

Under the product assumption $mbblue(rho)(mgreen(arrow(g)), mred(arrow(r)), t) = mbgreen(rho_G)(mgreen(arrow(g)), t) dot.c mbred(rho_R)(mred(arrow(r)), t)$, we can study the evolution of each marginal intensity separately. Moreover, if $mbgreen(rho_G)$ and $mbred(rho_R)$ each evolve according to independent PDEs, the product structure is preserved: the proposing and accepting landscapes evolve autonomously.

== PDE regimes on the intensity

Classic partial differential equations describe canonical modes of intensity evolution:

=== Diffusion

The heat equation
$ frac(partial mbblue(rho), partial t) = D nabla^2 mbblue(rho) $
models spreading or mixing. An initially concentrated intensity diffuses outward, representing diversification or loss of specificity. In the product case, if $mbgreen(rho_G)$ diffuses, entities become less specialized in their proposing behavior over time.

=== Advection

The transport equation
$ frac(partial mbblue(rho), partial t) = - arrow(v) dot nabla mbblue(rho) $
models directed drift. The intensity translates through the latent space at velocity $arrow(v)$, representing systematic change in interaction propensities. In ecological terms, this could model adaptation or environmental pressure shifting species' niches.

=== Reaction-diffusion

Combining local dynamics with spatial spreading:
$ frac(partial mbblue(rho), partial t) = D nabla^2 mbblue(rho) + f(mbblue(rho)) $
where $f(mbblue(rho))$ captures local growth, decay, or competition. This can produce pattern formation, traveling waves, or stable heterogeneous distributions.

== Induced dynamics on graph statistics

As $mbblue(rho)$ evolves, so do the expected graph properties. Under the product assumption, the expected number of nodes $bb(E)[N(t)] = mgreen(c_G)(t) dot.c mred(c_R)(t)$ and the expected edges $bb(E)[|E(t)|]$ become functions of time, determined by the evolving marginal intensities.

For instance, under pure diffusion with no-flux boundary conditions on the simplex, total mass is conserved: $mgreen(c_G)(t) = mgreen(c_G)(0)$. But the intensity-weighted means $mbgreen(mu_G)(t)$ and $mbred(mu_R)(t)$ may change, affecting expected edge counts even as expected node counts remain constant.

The inverse problem — inferring the PDE dynamics from a sequence of observed graphs $G_(t_1), G_(t_2), dots$ — combines the IDPG inference problem at each snapshot with dynamical estimation across time. This connects random graph theory to the broader literature on PDE inference from stochastic observations.

= An ecological motivating example<sec:foodwebs>

A _food web_ is an epitome example of an ecological network in which nodes represent species and edges represent consumption, predation, or, in brief, _who eats whom_. Usually, the set of species in a food web represents a certain ecological community or an ecosystem. And the edges are often determined by painstaking laborious field work by ecologists.

Yet, despite their fruitful application, it is quite clear that it is not a species eating another species: it's a certain individual of a species, say a cow, eating one or more individuals of another species. And, from a Darwinian point of a view, a species is a population (we might say a collection satisfying certain genealogical conditions) of individuals. Crucially, the individuals are not all identical, as various mechanisms bring some variance in the genetic (and phenotypical, and thus ecological) identity of individuals.

Hence, it is rather common in evolutionary sciences to describe the variability of individuals in a species with a certain probability distribution $mu$ in some space where the metric represents genetic similarity (and often the distributions tend to be multivariate Gaussians: most individuals will have a genome very similar to each other with few mutations, some will vary more, a few will be rather atypical, ...).

The environment, its biotic and abiotic components, and the phenotypic expression of an individual concur to determine the individual ecological role in an ecosystem (that is, the individual propensities to establish different, ecologically relevant, connections with other individuals from the same or other species.) This mapping of a genome to an ecological role corresponds to a mapping from the distribution $mu$, which takes values in the genetic space, to a distribution $cal(E)(mu)$ taking values in a theoretical space of ecological roles. In the case of a food web, where the ecological interactions are trophic relationships (who is a food resource for whom, who is a consumer of whom) we can ideally project $cal(E)(mu)$ into two subspaces $cal(E)_g(mu)$ and $cal(E)_r(mu)$ of _ecological role as a resource_ and _ecological role as a consumer_ (see @dalla2016exploring for such an analysis, although at the level of species).

An interaction between two individuals of different species will hence depend on the probability of those two individuals being "there", the propensity of one individual to eat the other, and the propensity of the latter to be eaten.

We can thus represent an ecological network, and in particular a food web, as an IDPG under the _edge-centric_ interpretation. The intensity $mbblue(rho)(mgreen(arrow(g)), mred(arrow(r)))$ captures the rate of potential trophic interactions between individuals with resource-role $mgreen(arrow(g))$ and consumer-role $mred(arrow(r))$. In the product case, $mbgreen(rho_G)$ and $mbred(rho_R)$ are given by (some transformation of) $cal(E)_g(mu)$ and $cal(E)_r(mu)$, the mapping of species genetic distributions into the ecological trophic-role spaces.

In this sense, a classic representation of a food web should be seen as a statistical summary of an edge-centric IDPG in which individual interaction sites are aggregated through an appropriate clustering procedure into species nodes. The necessity of considering food webs as probabilistic in nature has been recognized by @poisot2016structure.

Interestingly, ecological and evolutionary processes do really happen at the level of the underlying intensity functions, by the gradual movement of a species across the space, as fully recognized by various disciplines as population genetics or adaptive dynamics. The timescale of these evolutionary changes relative to the persistence of individual organisms determines whether the node-centric or edge-centric view is more appropriate for a given analysis.


#bibliography("bibliography.bib", style: "american-physics-society")


#pagebreak()

= Appendix: Derivations of expected edge counts<appendix:derivations>

We derive expected edge counts under the product assumption $mbblue(rho)(mgreen(arrow(g)), mred(arrow(r))) = mbgreen(rho_G)(mgreen(arrow(g))) dot.c mbred(rho_R)(mred(arrow(r)))$.

Recall the definitions:
- Marginal total intensities: $mgreen(c_G) = integral_mgreen(G) mbgreen(rho_G)(mgreen(arrow(g))) dif mgreen(arrow(g))$ and $mred(c_R) = integral_mred(R) mbred(rho_R)(mred(arrow(r))) dif mred(arrow(r))$
- Intensity-weighted mean positions: $mbgreen(mu_G) = integral_mgreen(G) mgreen(arrow(g)) mbgreen(rho_G)(mgreen(arrow(g))) dif mgreen(arrow(g))$ and $mbred(mu_R) = integral_mred(R) mred(arrow(r)) mbred(rho_R)(mred(arrow(r))) dif mred(arrow(r))$
- Normalized means: $mbgreen(tilde(mu)_G) = mbgreen(mu_G) \/ mgreen(c_G)$ and $mbred(tilde(mu)_R) = mbred(mu_R) \/ mred(c_R)$
- Expected number of nodes: $bb(E)[N] = mgreen(c_G) dot.c mred(c_R)$

== Node-centric derivation

The expected edge count is:
$ bb(E)[|E|] = integral_Omega integral_Omega (mgreen(arrow(g)) dot mred(arrow(r'))) lambda(dif mgreen(arrow(g)), dif mred(arrow(r))) lambda(dif mgreen(arrow(g')), dif mred(arrow(r'))) $

Under the product assumption:
$ bb(E)[|E|] = integral integral integral integral (mgreen(arrow(g)) dot mred(arrow(r'))) mbgreen(rho_G)(mgreen(arrow(g))) mbred(rho_R)(mred(arrow(r))) mbgreen(rho_G)(mgreen(arrow(g'))) mbred(rho_R)(mred(arrow(r'))) dif mgreen(arrow(g)) dif mred(arrow(r)) dif mgreen(arrow(g')) dif mred(arrow(r')) $

Since the dot product $mgreen(arrow(g)) dot mred(arrow(r')) = sum_k mgreen(g)_k mred(r')_k$ involves only $mgreen(arrow(g))$ and $mred(arrow(r'))$, the integral factorizes:
$ bb(E)[|E|] &= sum_k [integral mgreen(g)_k mbgreen(rho_G)(mgreen(arrow(g))) dif mgreen(arrow(g))] dot.c [integral mbred(rho_R)(mred(arrow(r))) dif mred(arrow(r))] dot.c [integral mbgreen(rho_G)(mgreen(arrow(g'))) dif mgreen(arrow(g'))] dot.c [integral mred(r')_k mbred(rho_R)(mred(arrow(r'))) dif mred(arrow(r'))] \
&= sum_k (mbgreen(mu_G))_k dot.c mred(c_R) dot.c mgreen(c_G) dot.c (mbred(mu_R))_k \
&= mgreen(c_G) dot.c mred(c_R) dot.c (mbgreen(mu_G) dot mbred(mu_R)) $

Rewriting in terms of normalized means:
$ bb(E)[|E|] = mgreen(c_G) dot.c mred(c_R) dot.c mgreen(c_G) dot.c mred(c_R) dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) = (bb(E)[N])^2 dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) $

== Edge-centric derivation

The edge intensity is $mbblue(rho_E)(mgreen(arrow(g)), mred(arrow(r))) = mbgreen(rho_G)(mgreen(arrow(g))) dot.c mbred(rho_R)(mred(arrow(r))) dot.c (mgreen(arrow(g)) dot mred(arrow(r)))$.

The expected edge count:
$ bb(E)[L] &= integral_mgreen(G) integral_mred(R) mbgreen(rho_G)(mgreen(arrow(g))) dot.c mbred(rho_R)(mred(arrow(r))) dot.c (mgreen(arrow(g)) dot mred(arrow(r))) dif mgreen(arrow(g)) dif mred(arrow(r)) \
&= sum_k [integral mgreen(g)_k mbgreen(rho_G)(mgreen(arrow(g))) dif mgreen(arrow(g))] dot.c [integral mred(r)_k mbred(rho_R)(mred(arrow(r))) dif mred(arrow(r))] \
&= sum_k (mbgreen(mu_G))_k dot.c (mbred(mu_R))_k = mbgreen(mu_G) dot mbred(mu_R) $

Rewriting:
$ bb(E)[L] = mgreen(c_G) dot.c mred(c_R) dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) = bb(E)[N] dot.c (mbgreen(tilde(mu)_G) dot mbred(tilde(mu)_R)) $
where we use $bb(E)[N] = mgreen(c_G) dot.c mred(c_R)$ as the common scale parameter.
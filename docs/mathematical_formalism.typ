#set document(title: "IDPG Mathematical Formalism", author: "IDPG Project")
#set page(margin: 2.5cm)
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")

#align(center)[
  #text(size: 20pt, weight: "bold")[Intensity Dot Product Graphs]
  #v(0.5em)
  #text(size: 14pt)[Mathematical Framework and Example Configurations]
]

#v(1em)

= Core Framework

== Latent Space

The latent space is the non-negative unit ball in $RR^d$:
$ B^d_+ = {bold(x) in RR^d : bold(x) >= 0, norm(bold(x)) <= 1} $

For $d = 2$, this is the first quadrant of the unit disk. For $d = 4$, it is the positive orthant of the 4-ball.

Each interaction site has two positions:
- $bold(g) in B^d_+$: source/giving position (propensity to send edges)
- $bold(r) in B^d_+$: receiver/target position (propensity to receive edges)

== Intensity Function

The intensity function $rho: B^d_+ times B^d_+ -> RR_+$ defines the density of interaction sites.

=== Product Intensity (ProductOfMixtures)

For general IDPG models, we use factorized intensities:
$ rho(bold(g), bold(r)) = rho_G (bold(g)) dot rho_R (bold(r)) $

Each marginal is a *BdPlusMixture* with $K$ components:
$ rho_G (bold(g)) = c_G sum_(k=1)^K w_k dot f_k (bold(g); bold(mu)_k, kappa_k) $

where:
- $c_G$: total intensity (scale)
- $w_k$: mixture weights ($sum_k w_k = 1$)
- $bold(mu)_k in B^d_+$: component mean positions
- $kappa_k > 0$: concentration parameters
- $f_k$: radial-angular density (von Mises-Fisher on direction, Beta on radius)

In this formulation, $bold(g)$ and $bold(r)$ are sampled *independently* from their respective marginals.

=== Mixture of Products (MixtureOfProductIntensities)

For ecological models with distinct species, we use a *sum of species-specific products*:
$ rho(bold(g), bold(r)) = sum_(m=1)^M rho_(G,m) (bold(g)) dot rho_(R,m) (bold(r)) $

Each species $m$ has:
- $rho_(G,m)$: intensity over G-space (resource/prey role)
- $rho_(R,m)$: intensity over R-space (consumer/predator role)

*Key differences from ProductIntensity:*

#table(
  columns: 2,
  stroke: 0.5pt,
  table.header([ProductOfMixtures], [MixtureOfProducts]),
  [$[sum_m rho_(G,m)] times [sum_m rho_(R,m)]$], [$sum_m [rho_(G,m) times rho_(R,m)]$],
  [$bold(g)$ and $bold(r)$ independent], [$bold(g)$ and $bold(r)$ coupled by species],
  [Cross-species mixing in $(bold(g), bold(r))$], [No cross-species mixing],
)

*Species abundance emerges from intensity:* The relative abundance of species $m$ is determined by:
$ gamma_m = c_(G,m) dot c_(R,m) = integral rho_(G,m) d bold(g) dot integral rho_(R,m) d bold(r) $

No extra abundance parameter is needed.

*Sampling from MixtureOfProducts:*
1. Total intensity: $C = sum_m gamma_m$
2. Sample $N tilde "Poisson"(C)$
3. For each site:
   - Select species $m$ with probability $gamma_m \/ C$
   - Sample $bold(g)$ from $rho_(G,m)$ (normalized)
   - Sample $bold(r)$ from $rho_(R,m)$ (normalized)
   - The site carries species identity $m$

== Poisson Point Process Sampling

Interaction sites are sampled from a Poisson Point Process (PPP) with intensity $rho$:
$ {(bold(g)_i, bold(r)_i)}_(i=1)^N tilde "PPP"(rho) $

The number of sites $N$ is Poisson distributed:
$ NN[N] = c_G dot c_R = integral_(B^d_+) rho_G (bold(g)) d bold(g) dot integral_(B^d_+) rho_R (bold(r)) d bold(r) $

== Connection Mechanism

Edge probability is the dot product:
$ PP(i -> j) = bold(g)_i dot bold(r)_j $

This creates directed graphs where:
- High $norm(bold(g))$ means strong sender
- High $norm(bold(r))$ means strong receiver
- Aligned directions increase connection probability

== Two Interpretations

*Node-Centric* (long-lived entities):
- Sampled sites become nodes
- Edges form independently with $PP(i -> j) = bold(g)_i dot bold(r)_j$
$ NN[|E|] = NN[N]^2 dot (tilde(bold(mu))_G dot tilde(bold(mu))_R) $

*Edge-Centric* (ephemeral interactions):
- Sample source site $S = (bold(g)_S, bold(r)_S)$ and target site $T = (bold(g)_T, bold(r)_T)$
- Connection probability: $PP(S -> T) = bold(g)_S dot bold(r)_T$
- All 4 coordinates are preserved per edge:
  - $bold(g)_S$: source's sending propensity (used for connection)
  - $bold(r)_S$: source's receiving propensity (available for clustering)
  - $bold(g)_T$: target's sending propensity (available for clustering)
  - $bold(r)_T$: target's receiving propensity (used for connection)
$ NN[L] = NN[N] dot (tilde(bold(mu))_G dot tilde(bold(mu))_R) $

*Full Site Clustering:* To identify entities, cluster based on the full $(bold(g), bold(r))$ signature (a $2d$-dimensional vector) rather than just $bold(g)$ or $bold(r)$ alone.

where $tilde(bold(mu))_G = bold(mu)_G \/ c_G$ is the normalized mean.

#pagebreak()

= Example Configurations

== Basic IDPG (2D)

*File:* `examples/basic_idpg.jl`

Source intensity $rho_G$ (2-component mixture):
#table(
  columns: 4,
  table.header([Component], [Weight], [Mean $bold(mu)$], [Concentration $kappa$]),
  [1], [0.6], [$(0.8, 0.2)$], [10.0],
  [2], [0.4], [$(0.2, 0.8)$], [10.0],
)
Total intensity: $c_G = 50.0$

Target intensity $rho_R$ (2-component mixture):
#table(
  columns: 4,
  table.header([Component], [Weight], [Mean $bold(mu)$], [Concentration $kappa$]),
  [1], [0.5], [$(0.3, 0.7)$], [10.0],
  [2], [0.5], [$(0.5, 0.5)$], [5.0],
)
Total intensity: $c_R = 50.0$

Expected values:
- $NN[N] = 50 times 50 = 2500$ interaction sites
- Average connection probability: $tilde(bold(mu))_G dot tilde(bold(mu))_R approx 0.37$

== IDPG vs Erdos-Renyi Comparison

*File:* `examples/idpg_vs_erdos_renyi.jl`

Three configurations tested:

*Configuration 1: Concentrated (orthogonal)*
$ rho_G: quad bold(mu) = (0.65, 0.15), quad kappa = 40, quad c = 80 $
$ rho_R: quad bold(mu) = (0.15, 0.65), quad kappa = 40, quad c = 80 $

*Configuration 2: Spread (aligned)*
$ rho_G: quad bold(mu) = (0.5, 0.5), quad kappa = 5, quad c = 40 $
$ rho_R: quad bold(mu) = (0.5, 0.5), quad kappa = 5, quad c = 40 $

*Configuration 3: Bimodal G*
$ rho_G: quad bold(mu)_1 = (0.8, 0.1), bold(mu)_2 = (0.1, 0.8), quad kappa = 30, quad c = 60 $
$ rho_R: quad bold(mu) = (0.5, 0.5), quad kappa = 10, quad c = 60 $

Key finding: IDPG degree variance/mean ratio is 2-20x higher than Erdos-Renyi.

#pagebreak()

= Ecological Food Web Models

Ecological food webs use *MixtureOfProductIntensities* because each species has a coupled $(bold(g), bold(r))$ niche:
$ rho(bold(g), bold(r)) = sum_(m=1)^M rho_(G,m) (bold(g)) dot rho_(R,m) (bold(r)) $

When sampling, both $bold(g)$ and $bold(r)$ come from the *same species* - no cross-species mixing.

== 2D Food Web

*File:* `examples/ecological_example.jl`

Using hyperspherical coordinates $(r, phi)$ where $r$ is radius and $phi in [0, pi/2]$ is angle:
$ bold(x) = r dot (cos phi, sin phi) $

Three species (trophic guilds), each with coupled G and R niches:
#table(
  columns: 5,
  table.header([Species $m$], [$(r, phi)_G$], [$(r, phi)_R$], [G-role], [R-role]),
  [Producers], [$(0.85, pi/10)$], [$(0.02, pi/6)$], [Visible resource], [Non-consumer],
  [Herbivores], [$(0.85, pi/4)$], [$(0.85, pi/5)$], [Mid-chain resource], [Consumes producers],
  [Predators], [$(0.20, pi/3)$], [$(0.85, pi/4)$], [Rarely eaten], [Consumes herbivores],
)

Each species $m$ has intensity $rho_m (bold(g), bold(r)) = rho_(G,m)(bold(g)) dot rho_(R,m)(bold(r))$ with concentration $kappa = 50$.

== 4D Food Web

*File:* `examples/ecological_4d_example.jl`

Using 4D hyperspherical coordinates $(r, phi_1, phi_2, phi_3)$:
$ bold(x) = r dot (cos phi_1, sin phi_1 cos phi_2, sin phi_1 sin phi_2 cos phi_3, sin phi_1 sin phi_2 sin phi_3) $

Angle constants for trophic separation:
$ phi_"small" = pi/60 approx 3 degree, quad phi_"large" = 29pi/60 approx 87 degree $

Five species with coupled $(bold(g), bold(r))$ niches:

#table(
  columns: 5,
  table.header([Species $m$], [$r_G$], [$r_R$], [G-niche (resource role)], [R-niche (consumer role)]),
  [Producers], [0.95], [0.02], [Dim 1: highly visible], [Nearly zero consumption],
  [Small Herb.], [0.90], [0.95], [Dim 2: herbivore prey], [Dim 1: targets producers],
  [Large Herb.], [0.85], [0.90], [Dim 2: herbivore prey], [Dim 1: targets producers],
  [Small Pred.], [0.75], [0.90], [Dim 3: predator prey], [Dim 2: targets herbivores],
  [Apex Pred.], [0.20], [0.95], [Dim 4: rarely eaten], [Dim 3: targets small pred.],
)

The intensity for species $m$ is:
$ rho_m (bold(g), bold(r)) = rho_(G,m)(bold(g); bold(mu)_(G,m), kappa) dot rho_(R,m)(bold(r); bold(mu)_(R,m), kappa) $

Parameters: $kappa = 30$ for all species.

*Species abundance* is determined by $gamma_m = c_(G,m) dot c_(R,m)$, with no separate abundance parameter.

#pagebreak()

= PDE Evolution

== Diffusion on $B^d_+$

*File:* `examples/diffusion_example.jl`

The intensity evolves according to the heat equation:
$ (diff rho) / (diff t) = D nabla^2 rho $

with absorbing boundary conditions at $norm(bold(x)) = 1$.

*Parameters:*
- Initial: Gaussian at $bold(mu) = (0.8, 0.2)$ with $kappa = 15$, $c = 50$
- Diffusion coefficient: $D = 0.1$
- Time step: $Delta t = 0.0001$
- Final time: $t_"final" = 2.0$

Effect: Intensity spreads and total mass decreases (absorbed at boundary).

== Advection on $B^d_+$

Intensity transported by velocity field $bold(v)$:
$ (diff rho) / (diff t) + nabla dot (rho bold(v)) = 0 $

For constant velocity, this shifts the distribution.

== Temporal Food Web Dynamics

*File:* `examples/temporal_foodweb.jl`

Four regimes compared:

*1. Static:* No evolution (baseline)

*2. Diffusion:*
$ D = 0.08, quad Delta t = 0.002, quad t_"final" = 1.0 $

*3. Advection:*
$ bold(v) = (0.4, 0.3, -0.1, 0.0), quad Delta t = 0.002, quad t_"final" = 1.0 $

*4. Pursuit-Evasion (coupled dynamics):*

Resources flee consumers:
$ bold(v)_G (bold(x)) = v_"pursuit" dot (bold(x) - tilde(bold(mu))_R) / norm(bold(x) - tilde(bold(mu))_R) + v_"center" dot (bold(x)_"center" - bold(x)) $

Consumers chase resources:
$ bold(v)_R (bold(x)) = v_"pursuit" dot (tilde(bold(mu))_G - bold(x)) / norm(tilde(bold(mu))_G - bold(x)) + v_"center" dot (bold(x)_"center" - bold(x)) $

where:
- $v_"pursuit" = 0.5$ (chase/flee speed)
- $v_"center" = 0.15$ (centering force)
- $bold(x)_"center" = (0.25, 0.25, 0.25, 0.25)$ (center of $B^4_+$)

Grid: resolution 12 ($12^4 = 20736$ points, filtered to $approx 700$ inside $B^4_+$)

#pagebreak()

= Key Formulas Summary

== Expected Quantities (Product Intensity)

#align(center)[
#table(
  columns: 3,
  stroke: 0.5pt,
  table.header([Quantity], [Formula], [Description]),
  [$c_G, c_R$], [$integral rho_G d bold(g)$, $integral rho_R d bold(r)$], [Marginal intensities],
  [$bold(mu)_G, bold(mu)_R$], [$integral bold(g) rho_G d bold(g)$, $integral bold(r) rho_R d bold(r)$], [Intensity-weighted means],
  [$tilde(bold(mu))_G, tilde(bold(mu))_R$], [$bold(mu)_G / c_G$, $bold(mu)_R / c_R$], [Normalized means],
  [$NN[N]$], [$c_G dot c_R$], [Expected sites],
  [avg. conn. prob.], [$tilde(bold(mu))_G dot tilde(bold(mu))_R$], [Mean $PP(i -> j)$],
  [$NN[|E|]$], [$NN[N]^2 dot (tilde(bold(mu))_G dot tilde(bold(mu))_R)$], [Node-centric edges],
  [$NN[L]$], [$NN[N] dot (tilde(bold(mu))_G dot tilde(bold(mu))_R)$], [Edge-centric edges],
)
]

== Expected Quantities (Mixture of Products)

For $rho(bold(g), bold(r)) = sum_m rho_(G,m)(bold(g)) dot rho_(R,m)(bold(r))$:

#align(center)[
#table(
  columns: 3,
  stroke: 0.5pt,
  table.header([Quantity], [Formula], [Description]),
  [$c_(G,m), c_(R,m)$], [$integral rho_(G,m) d bold(g)$, $integral rho_(R,m) d bold(r)$], [Per-species marginal intensities],
  [$gamma_m$], [$c_(G,m) dot c_(R,m)$], [Species $m$ total intensity],
  [$C$], [$sum_m gamma_m$], [Total intensity (all species)],
  [$NN[N]$], [$C = sum_m gamma_m$], [Expected sites],
  [$PP("species" = m)$], [$gamma_m \/ C$], [Probability site is from species $m$],
)
]

*Key difference:* In MixtureOfProducts, sites carry species identity. The $(bold(g), bold(r))$ pair comes from the same species, enabling species-level analysis.

== Relationship Between Interpretations

$ frac(NN[|E|], NN[L]) = NN[N] $

The node-centric interpretation produces $NN[N]$ times more edges because each of the $N$ nodes can connect to any other node.

== Mesoscale Properties

- *Reciprocity* $approx tilde(bold(mu))_G dot tilde(bold(mu))_R$ (edges are independent)
- *In/Out degree correlation* $approx 0$ (for product intensity with independent $G$ and $R$)
- *Degree variance/mean* $>>$ 1 (distinguishes IDPG from Erdos-Renyi)
- For MixtureOfProducts: *In/Out correlation* can be non-zero due to species-level coupling

== Volume of $B^d_+$

$ "Vol"(B^d_+) = frac(pi^(d\/2), 2^d dot Gamma(d/2 + 1)) $

For $d = 2$: $pi/4 approx 0.785$. For $d = 4$: $pi^2/32 approx 0.308$.

Higher dimensions require larger intensity scales to sample similar numbers of sites.

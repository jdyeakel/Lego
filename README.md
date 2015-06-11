## Lego Anime model

Community assembly and evolution via discrete interaction rules


This is a work in progress.

```S
A: avoid
N: need
I: ignore
M: make
E: eat
```

Initial matrix seeding rules

Sp A | Sp B
--- | ---
E | {I, M, E}
M | N
I | {A, N, I, E}
N | {N, I, M, E}
A | {A, I}

Biological mapping

Interaction | Biological meaning
--- | ---
{N,E} | Facultative mutualism
{N.N} | Obligate mutualism
{N,I} | Commensalism
{I,A} | Asymmetric unstable comp.
{A,A} | Symmetric unstable comp.
{I,E} | Asymmetric Predation
{M,N} | Engineering
{I,I} | Stable competition, coexistence
{E,E} | Symmetric predation
{N,A} | DNE
{E,A} | DNE
{M,M} | DNE
{M,A} | DNE
{M,I} | DNE
{M,E} | DNE



on diagonal CAN determine active, but not passive


RULES:
1. if a -M-> b, then b -I-> everything, except b -N-> a

2. if similarity(a,b) > X; then a -A-> b

3. active players must consume: number(E) > 0

4. row 1 is the sun, R1 = I; C1 = {E,I}

RULES FOR CULLING:
1. if nothing from a given set ?(primary producers) is present; then everything else is culled














##Old version


The goal of the `Lego` project is to provide a fundamental model for the evolution of biodiversity and organismal niche spaces. The `Lego` project is currently built around a mass-balance equation
```S
A + B = C
```
where `A` represents a species, `B` is the set of coproducts produced by the species, and `C` is the set of resources required by `A`. 
Initially (i.e. for the first eden species), the set of `Global Resources` includes `R = {C, B, A}`.

`Speciation` occurs when `A` => `A'`, with coproducts `B'` and `C'` based on some mutation rate. 
The new resources `C'` are modified from `C`, and partly drawn from the `Global Resources`.
The coproducts `B'` are modified from `B`, and partly drawn from the `Global Resources`.
There is some probability of creating a new coproduct.
The `Global Resources` are then updated to include `R = {C, C', B, B', A, A'}`.
The similarity of `A` and `A'` is measured as a function of the similarity between `C` and `C'`.

From these basic reactions, the following ecological relationships can be observed:
```
Competition: (A vs. A') determined as a function of Sim(C,C').
Facilitation: When B is in the C' of A', then the persistence of A facilitates the persistence of A'.
Predation: When the C' of A' includes A.
Engineering: When the B' of A' has high overlap with the resources of many species.
```

The PseudoCode:


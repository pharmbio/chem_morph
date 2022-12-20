<style>
#box {
  background-color: #fffad8;
  border-radius: 3pt;
  border: 2px solid #a2a2a2;
  padding: 10px;
  margin-top: 10px;
  margin-bottom: 10px;
  margin-right: 5px;
  margin-left: 5px;
}
</style>
------------------------------------------------------------------------

**ANALYSIS REPORT**  

A survey on different approaches to calculate distance in cell profiling fingerprint space and chemical space and their relationship


------------------------------------------------------------------------

![#f03c15](https://placehold.it/15/f03c15/000000?text=+)

**Pharmaceutical Bioinformatics Research Group**  

by Nima Chamyani  


# Introduction

Every now and then, we encounter a situation in which we wish to
measure, determine or quantify relationships between some variables by
observing their various characteristics, that may or may not be
correlated. Numerous mathematical methods can be used to find these
relationships, but they only work when they are being used for the very
specific question they tend to answer. It is also possible to
misinterpret some similar concepts, such as similarity and correlation,
where they represent different things. To make thing more complicated,
it can be challenging to compare the similarities that have been
obtained on different spaces with completely different descriptors.

For long the assumption of chemical similarity may result in similar
activities has been a foundation for developing different similarity
algorithms which can use topological, 3D conformation and/or physical
properties. In contrast, in cell profiling we just have a set of
independent compound descriptors to work with. To determine the
similarity between two chemical spaces, there are a number of widely
used methods. They will encompass different structural properties and I
will go through them in the method section. Furthermore, I will explore
different mathematical distance metrics and their main purpose to
perform on cell profile features too.

# Aim of the analysis

In this analysis, I attempted to find the best way to calculate pairwise
similarity of 931 compounds (SPECS-935 dataset) both in chemical
structure space and in cell profiling finger print space, then compare
these two similarities to each other for evaluation of any relationship.
For this I will explore different distance measurement technique and I
will elaborate on my choice. At the end I will show the result and
discuss about the possible interpretation.

# Methods

## Distance metrics

It is meaningless to talk about \"close\" and \"far\" in any space, when
there is no measurable distance. In order to establish these notions
over a set of abstract mathematical objects, we must be capable of
assessing how close each pair of these objects is to each other and
interpret it as similarity. There are several different similarity
functions that are used for measuring the distance between two vectors,
numbers, or pairs of numbers. However, a proper distance measure should
have a few properties to be metric. To be precise, a distance metric is
a function with the following properties:

1.  If the distance of two objects is zero, then they are the same, and
    vice versa $$d(x,y) = 0 \Rightarrow x = y$$

2.  A proper distance metric is symmetric $$d(x,y)=d(y,x)$$

3.  A proper distance metric satisfies triangular inequality (meaning
    for any triangle connecting three point, the sum of the lengths of
    any two sides must be greater than or equal to the length of the
    remaining side) $$d(x,y) \le d(x,z)+d(z,y)$$

Now if we center the distance between 0 and 1 then the similarity (S)
can be define as $$S(x,y) = 1 - d(x,y)$$

### Correlation and distance correlation

***Correlation***: Correlation is a statistical measure that expresses
the extent to which two variables are linearly related. Since
correlation is looking for the trend on the samples, it is not a proper
distance metrics. It cannot follow the triangular inequality meaning
except for the extreme cases, it's not possible to find the correlation
of two random variables given their correlations with a third one. So
it's not possible to use correlation to measure the distances and
similarity. For instance, the correlation between $(1, 2, 3)$ and
$(10, 20, 30)$ or any other set with the same trend and different
magnitude like $(100, 200, 300)$ is the same and equal to 1. It cannot
be interpreted that they are the same point in the space or close to
each other. But they act in a similar manner because they all seem to
point in the same direction.

There are several type of correlation algorithms. Pearson correlation is
the most widely used correlation statistic to measure the degree of the
relationship between linearly related variables. It assumes that
variables are normally distributed, have linearity between them and data
is equally distributed around the regression line (homoscedasticity).
Another type of correlation is Kendall rank correlation which is a
non-parametric test that measures the strength of dependence between two
variables. While linear relationships mean two variables move together
at a constant rate (i.g. two parallel and straight line), monotonic
relationships measure how likely it is for two variables to move in the
same direction, but not necessarily at a constant rate (not a linear
path). Spearman's is incredibly similar to Kendall's. It is a
non-parametric test that measures a monotonic relationship using ranked
data. While it can often be used interchangeably with Kendall's,
Kendall's is more robust and generally the preferred method of the
two(Table [1](#corr){reference-type="ref" reference="corr"}). The
general formula for correlation is:

$$corr(x,y) = \frac{cov(X, Y)}{ \sqrt{ \sigma_{X} \times \sigma_{Y} }}$$

where $cov(X, Y)$ denotes the covariance of $X$ and $Y$ and $\sigma$ is
the variance of the variable.

::: {#corr}

| Correlation | Type of relation | Type of measurement                        | Type of distribution |
|-------------|------------------|--------------------------------------------|----------------------|
| Pearson     | Linear           | Quantitative (interval or ratio) variables | Normal               |
| Spearman    | Non-linear       | Ordinal, interval or ratio variables       | Any                  |
| Kendall     | Non-linear       | Ordinal, interval or ratio variables       | Any                  |

  Table 1 : Most widely used correlation techniques

:::

***Distance correlation***: The correlation value can range from -1 to
1, which indicates perfect negative or positive correlation while zero
represents no correlation at all. Also, in correlation the converse
implication is not true, meaning the result of reversing its two
constituent statements is not valid. In another term, zero correlation
between two variable does not necessarily imply their independence. This
drawback lead to developing a new methods that is called distance
correlation. So beside working with non-linear type of relationship, a
distance correlation of zero indicates that there is in fact no
dependence between the two variables.

$$dCorr^2(x,y) = \frac{dCov^2(X, Y)}{ \sqrt{ d\sigma^2_{X} \times d\sigma^2_{Y} }}$$

The distance correlation produces non-negative value between 0 and 1.

### Euclidean distance

The Euclidean distance between two points measures the length of the
shortest segment connecting them. It is the most obvious way of
representing distance between two points that can be calculated by the
Pythagorean theorem.
$$d(x, y)=\sqrt{\sum_{i=1}^k\left(x_i-y_i\right)^2}$$

The Euclidean distance is most commonly used to calculate the distance
between two rows of numerical data, such as floating-point or integer
data. It is often done after normalizing or standardizing numeric values
otherwise, the distance measure will be dominated by large values.

::: {#box}
NOTE: In order to speed up distance calculations, it is common to remove
the square root operation when performing thousands or millions of
calculations. After this modification, the scores will still have the
same relative proportions and can still be used effectively within a
machine learning algorithm.
:::

### Manhattan distance

In Manhattan distance, the distance between two points is measured by
the absolute sum of their coordinate differences (also referred to as
rectilinear distance, L1 distance, taxicab geometry or city block
distance).

$$d_m(\mathbf{x}, \mathbf{y})=\|\mathbf{x}-\mathbf{y}\|_m=\sum_{i=m}^n\left|x_i-y_i\right|$$

It might make sense to calculate Manhattan distance instead of Euclidean
distance for two vectors in an integer feature space where vectors
describe objects on a uniform grid. Hence the taxicab name for the
measure: the shortest path a taxicab would take between city blocks
(coordinates on the grid).

### Minkowski distance

Minkowski distance calculates the distance between two real-valued
vectors. It is generalized from of Euclidean and Manhattan methods and
by adding \"p\" as a parameter, it allows different distance measures to
be calculated.

$$D(X, Y)=\left(\sum_{i=1}^n\left|x_i-y_i\right|^p\right)^{1 / p}$$

So if p is 2, Minkowski distance is the same as the Euclidean distance
and when p is 1, it's the same as the Manhattan distance, however,
intermediate values provide a controlled balance between the two
measures. This is particularly useful in machine learning algorithms
that use distance measures, since it gives control over the type of
distance measure to be used for real-valued vectors.

### Cosine similarity

The cosine similarity computes the similarity between two samples that
have been obtained from the same or different distributions. Samples are
viewed as vectors in an inner product space, and the cosine similarity
is defined as the cosine of the angle between them, that is, the dot
product of the vectors divided by the product of their lengths.

$$S_c(x,y) = cos(\theta) = \frac{\sum_i x_i y_i}{ \sqrt{ \sum_i x_i^2} \sqrt{ \sum_i y_i^2 } } 
= \frac{ \langle x,y \rangle }{ ||x||\ ||y|| }$$

As a result, cosine similarity is used when the magnitude between two
vectors is not relevant and it is only their orientation that will
determine their similarity. Hence, two vectors with the same orientation
have a cosine similarity of 1, two vectors at 90° have a similarity of
0, and two vectors diametrically opposed have a similarity of -1. Thus,
cosine similarity resembles Pearson correlation very closely. In a
similar way to Pearson correlation, the cosine distance cannot be
considered as a proper distance metric due to the fact that it does not
exhibit the triangle inequality property. It is therefore very common to
run into confusion when comparing the cosine similarity and the Pearson
correlation. However, correlation is just acting similar to the cosine
similarity where vectors are being centred ($\bar{x} = \bar{y} = 0$) and
unlike the cosine, the correlation is invariant not only to scale but
also to the shift (location changes) of x and y.

::: {#box}
**Cosine Similarity vs Pearson Correlation**\
We know that the cosine similarity between two vectors $a$ and $b$ is
just the angle between them
$$\cos\theta = \frac{a\cdot b}{\lVert{a}\rVert \, \lVert{b}\rVert}$$

And for $a$ vector $x$ the \"$z$-score\" vector would typically be
defined as $$z=\frac{x-\bar{x}}{s_x}$$ where
$\bar{x}=\frac{1}{n}\sum_ix_i$ and $s_x^2=\overline{(x-\bar{x})^2}$ are
the mean and variance of $x$. So $z$ has mean 0 and standard deviation
1, i.e. $z_x$ is the standardized version of $x$. For two vectors $x$
and $y$, their correlation coefficient would be
$$\rho_{x,y}=\overline{(z_xz_y)}$$

Now if the vector a has zero mean, then its variance will be
$s_a^2=\frac{1}{n}\lVert{a}\rVert^2$, so its unit vector and $z$-score
will be related by
$$\hat{a}=\frac{a}{\lVert{a}\rVert}=\frac{z_a}{\sqrt n}$$

So if the vectors a and b are centered (i.e. have zero means), then
their cosine similarity will be the same as their correlation
coefficient.
:::

### Triangle Area Similarity -- Sector Area Similarity (TS-SS)

**Problems with Euclidean base metrics**: Two vectors with no common
component values may have a smaller distance than similar vectors
containing the same component values. For instance, in A(1, 0, 0), B(0,
-2, -2), C(10, 0, 0) if we calculate the distance we would have d(A , C)
= 9 units and d(A , B) = 3 units. So, according to Euclidean metric, A
and B are closer or more similar than A and C. Despite this, A and C are
in the same direction and only their magnitude differs. So trends in
vectors cannot be explained by Euclidean metric systems.

**Problems with Cosine similarity**: As it mentioned before, Cosine
similarity is invariant to scale. As a result, it does not take
magnitude into account and only focuses on orientation. As an example,
all three vectors of A(1, 2, 3), B(10, 20, 30) and C(100, 200, 300) will
be considered to be the same. Despite the fact that vector A and vector
B are closer to each other than any other combination of A, B, and C,
cosine similarity cannot further distinguish them.

Triangle Area Similarity -- Sector Area Similarity (TS-SS) was developed
to address this problem. An algorithm that can combine both the
direction and magnitude of vector in similarity check.

::: {#box}
NOTE: I found this method in an article originally falls under Natural
Language Processing (NLP) but the logic of the metric seems to work very
well for any vector similarity check.\
1- Euclidean distance:
$$ED(x, y)=\sqrt{\sum_{i=1}^k\left(x_i-y_i\right)^2}$$

2- Triangle's Area Similarity (TS):
$$\operatorname{TS}(x, y)=\frac{|x| \cdot|y| \cdot \sin \left(\theta^{\prime}\right)}{2}$$
3- The magnitude difference between two vectors:
$$\operatorname{MD}(x, y)=\left|\sqrt{\sum_{i=1}^k x_i^2}-\sqrt{\sum_{i=1}^k y_i^2}\right|$$
4- Sector's Area Similarity (SS):
$$\mathrm{SS}(x, y)=\pi \cdot(\mathrm{ED}(x, y)+\operatorname{MD}(x, y))^2 \cdot\left(\frac{\theta^{\prime}}{360}\right)$$

5- Finally by multiplying them together, we combine TS and SS:

$$d_{TS-SS} = TS \times SS$$
:::

### Jaccard-Needham similarity

In 1884, Grove Karl Gilbert developed the Jaccard similarity, which was
later independently developed by Paul Jaccard and Tanimoto. Their work
are identical in generally taking the ratio of intersection over union
of the sets. Thus, it measures the similarity between finite sample sets
by dividing the portion that they have in common minus the part that is
different.

$$J(X, Y)=\frac{|X \cap Y|}{|X \cup Y|}=\frac{|X \cap Y|}{|X|+|Y|-|X \cap Y|}$$

It is is widely used in computer science where binary or binarized data
are present. In confusion matrices employed for binary classification,
the Jaccard index can be framed in the following formula:

$${\displaystyle {\text{Jaccard index}}={\frac {TP}{TP+FP+FN}}}$$

where TP are the true positives, FP the false positives and FN the false
negatives

::: {#box}
NOTE: Various forms of functions described as Tanimoto similarity and
Tanimoto distance occur in the literature and on the Internet. Most of
these are synonyms for Jaccard similarity and Jaccard distance, but some
are mathematically different. The similarity ratio is equivalent to
Jaccard similarity, but the distance function is not the same as Jaccard
distance.
:::

### Sørensen--Dice similarity

Sørensen--Dice similarity, F1 score, Czekanowski's binary
(non-quantitative) index or Zijdenbos similarity index all referring to
an index equals twice the intersection (number of elements common to
both sets) divided by the sum of the elements in each set.

$${\displaystyle DSC(X, Y)={\frac {2|X\cap Y|}{|X|+|Y|}}}$$

When applied to Boolean data, using the definition of true positive
(TP), false positive (FP), and false negative (FN), it can be written
as:

$${\displaystyle DSC={\frac {2TP}{2TP+FP+FN}}}$$

::: {#box}
NOTE: Sørensen--Dice similarity is not very different in form the
Jaccard similarity. Both are equivalent in the sense that given a value
for the Sørensen--Dice coefficient ${\displaystyle S}$, one can
calculate the respective Jaccard index value ${\displaystyle J}$ and
vice versa, using the equations ${\displaystyle J=S/(2-S)}$ and
${\displaystyle S=2J/(1+J)}$.

Since the Sørensen--Dice coefficient does not satisfy the triangle
inequality, it can be considered a semimetric version of the Jaccard
index. The function ranges between zero and one, like Jaccard. Unlike
Jaccard, the corresponding difference function

$${\displaystyle d=1-{\frac {2|X\cap Y|}{|X|+|Y|}}}$$

is not a proper distance metric as it does not satisfy the triangle
inequality.
:::

## Chemical Similarity

### Constitutional (topological) similarity

The 2D molecular structure are the bases to assess topological
similarity. A number of approaches have been proposed to extract
information from molecular structures. However, these structural models
cannot distinguish conformers due to their inability to interpret 3D
molecular structures.

**Classic Topological Descriptors**: The most basic representation of
topology is the count of individual atoms, bonds, rings, pharmacophore
points, and the degree of connectivity between the atoms. Using these
two-dimensional fragment descriptors (atom-centered, bond-centered,
ring-centered fragments) the similarity of two compound will be
assessed.

**Molecular Fingerprints and Molecular Holograms**: In general, a
molecule's fingerprint can be interpreted as a system of encoding a
molecule's structure. The most common fingerprint is a sequence of
binary digits (bits) that indicate whether a molecule contains certain
substructures. So, by using bit-strings (fingerprints), two-dimensional
substructures can be encoded to a vector. Different distance algorithms
(such as those discussed in the distance metric section) can be applied
to these binary fingerprints in order to determine their distances.

Given two binary vectors, $X$ and $Y$, each with n binary attributes,
each attribute of $X$ and $Y$ can either be 0 or 1. The total number of
each combination of attributes for both $X$ and $Y$ are specified as
follows:

-   $a$ represents the total number of attributes where $X$ and $Y$ both
    have a value of 1.

-   $b$ represents the total number of attributes where the attribute of
    $X$ is 1 and the attribute of $Y$ is 0.

-   $c$ represents the total number of attributes where the attribute of
    $X$ is 0 and the attribute of $Y$ is 1.

-   $d$ represents the total number of attributes where $X$ and $Y$ both
    have a value of 0.

-   $n = a + b + c + d$

Different similarity algorithms is listed in Table
[2](#alg){reference-type="ref" reference="alg"} based on these
parameters.

::: {#alg}

| **Algorithm**              | **Formula**                                                                 |
|--------------------------------|-----------------------------------------------------------------------------------|
| Jaccard                        | $a / (a+b+c)$                                                                     |
| Russel - Rao                   | $a/n$                                                                             |
| Rogers - Tanimoto              | $(a+d)/(a+2 \times (b+c)+d)$                                                      |
| Kulczynski \#1                 | $a/(b+c)$                                                                         |
| Kulczynski \#2                 | $0.5 \times (a/(a+b)+a/(a+c))$                                                    |
| Dice                           | $2 \times a/(2 \times a+b+c)$                                                     |
| Pearson's Phi coefficient      | $((a \times d)-(c \times b))/\sqrt{(a+c) \times (c+d) \times (a+b) \times (b+d)}$ |
| Baroni-Urbani/Buser            | $(a+\sqrt{a \times d})/(a+b+c+\sqrt{a \times d})$                                 |
| Braun-Blanquet                 | if $(a+b)> (a+c)$  then $S = a/(a+b)$ else $S = a/(a+c)$                          |
| Simpson similarity coefficient | if $(a+b) < (a+c)$ then $S = a/(a+b)$ else $S = a/(a+c)$                          |
| Michael                        | $4 \times (a \times d-b \times c)/((a+d) \times (a+d)+(b+c) \times (b+c))$        |
| Sokal and Sneath \#1           | $a/(a+2 \times (b+c))$                                                            |
| SokalSneath \#2                | $0.25  \times  ( a/(a+b)+ a/(a+c)+ d/(b+d)+ d/(c+d) )$                            |
| SokalSneath \#3                | $a \times d/\sqrt{(a+b) \times (a+c) \times (d+b) \times (d+c)}$                  |
| Sokal and Sneath \#4           | $(a+d)/(b+c)$                                                                     |
| Sokal and Sneath \#5           | $2 \times (a+d)/(2 \times (a+d)+b+c)$                                             |
| Simple Matching                | $(a+d)/(a+b+c+d)=(a+d)/n$                                                         |
| Sneath - Sokal                 | $(a+d)/(a +0.5 \times (b+c)+d) $                                                  |
| Kocher - Wong                  | $a \times n/((a+b) \times (c+d))$                                                 |
| Ochiaï \#1                     | $a/\sqrt{(a+b) \times (a+c)}$                                                     |
| Ochiaï \#2                     | $a \times d/\sqrt{(a+b) \times (a+c) \times (d+b) \times (d+c)}$                  |
| Yule's Sigma                   | $(\sqrt{a \times d}-\sqrt{b \times c})/(\sqrt{a \times d}+\sqrt{b \times c})$     |
| Yule's Q                       | $(a \times d-b \times c)/(a \times d+b \times c)$                                 |
| McConnoughy                    | $(a \times a - b \times c) / \sqrt{(a+b) \times (a+c)}$                           |
| Phi Square                     | $(a \times d + b \times c)^2 / ((a+b) \times (a+c) \times (b+c) \times (b+d))$    |
| Dispersion                     | $(a \times d-b \times c)/(a+b+c+d)^2$                                             |

  Table 2 : Similarity algorithm for binary vectors.
  
:::

**Burden eigenvalue descriptors or BCUT**: BCUT descriptors are based on
an extension of Burden's approach[@burden1989molecular] for searching
large databases for chemical similarity[@pearlman2002novel]. A molecule
graph with hydrogen included is used to calculate the Burden matrix, and
condensing the information from those matrices to eigenvalues results in
a one-dimensional measure, which closely reflect the structure of a
molecule. In particular, BCUTs are widely used in similarity analyses.

### Configuration and conformation similarity

Some activities may require more consideration of conformational
flexibility of chemical compounds and their 3D descriptors such as shape
and volume not just topological information. It is possible, in these
cases, to model activity and evaluate similarity only with the molecular
conformations. This is widely used in virtual screening of the ligands
to see their binding affinity to proteins. The same ligand cavity in the
receptor can be filled with two molecules with very different
topological structures but very similar 3D conformations. Therefore,
it's safe to assume that conformation similarity checks should yield a
more realistic score than topological similarity tests for most of the
cases.

Configuration and conformation similarity checks can be performed using
a number of different approaches:

**Shape**: **Distance-based and Angle-based Descriptors** **Three
Dimensional and Field Similarity** **Molecular Multi-pole Moments**

### Physicochemical properties similarity

### Quantum-chemistry similarity

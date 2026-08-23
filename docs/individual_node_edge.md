---
kernelspec:
  name: python3
  display_name: python3
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.13'
    jupytext_version: 1.13.8
---

```{eval-rst}
.. currentmodule:: tstrait
```

(genetic_individual_node_edge_doc)=

# Individual, node, and edge genetic values

In line with the {ref}`tree sequence data model<tskit:sec_data_model>`,
this page describes relationships between
individual, node, and edge genetic values,
and underlying edge effects arising from mutations.

**Learning Objectives**

After this page, you will be able to:

- Understand how to obtain genetic value in tstrait for individuals, nodes, and edges
- Understand how causal-allele effects give rise to edge effects and
  the genetic values.

# Algorithm Overview

The tstrait algorithm for [individuals' genetic values](genetic_value_doc)
can also be used to compute node and edge genetic values,
and there is the related concept of edge effects;
all following the assumed [quantitative genetics model](phenotype_model).
These quantities are related as follows:

- *edge effect* is the sum of effects of causal alleles introduced on the edge by mutations.
- *edge "genetic" value* is its edge effect plus a sum of edge effects above the edge.
- *node "genetic" value* is a sum of "genetic" values of incoming edges into the node.
- *individual genetic value* is a sum of "genetic" values of nodes of the individual.

Above we have quoted the term "genetic" when referring to nodes and edges,
because that term is usually used for individuals.
Because individuals' genetic value is a sum over the contributions of its nodes,
the term "genetic value" applies also to nodes, and hence also to incoming edges.

# Example

We will demonstrate the concepts of edge effects and
edge, node, and individual genetic values with
a tiny example so we can follow the calculations.

```{code-cell}
import io
import tskit
import tstrait
import pandas as pd

individuals = io.StringIO(
    """\
    flags location parents
        0        0      -1
        0        0      -1
    """
)

nodes = io.StringIO(
    """\
    id is_sample time individual
     0         1    0          0
     1         1    0          0
     2         1    0          1
     3         1    0          1
     4         0    1         -1
     5         0    1         -1
     6         0    2         -1
     7         0    4         -1
    """
)

edges = io.StringIO(
    """\
    left right parent child
       0    10      4     0
       0     5      4     1
       5    10      5     1
       0    10      5     2
       0    10      6     3
       0    10      6     5
       0    10      7     4
       0    10      7     6
    """
)

sites = io.StringIO(
  """\
  position ancestral_state
         0               0
         1               0
         2               0
         3               0
         4               0
         5               0
         6               0
         7               0
         8               0
         9               0
  """
)

mutations = io.StringIO(
  """\
  site node time derived_state parent
     0    4  3.0             1     -1
     1    3  1.5             1     -1
     2    1  0.5             1     -1
     4    4  2.0             1     -1
     6    6  3.0             1     -1
     7    5  1.5             1     -1
     8    3  1.0             1     -1
     9    1  0.5             1     -1
  """
)

ts = tskit.load_text(individuals = individuals,
                     nodes = nodes,
                     edges = edges,
                     sites = sites,
                     mutations = mutations,
                     strict = False)
print(ts)

ts.genotype_matrix().T

data = {"site_id":       [  0,   1,   2,   4,   6,   7,   8,   9],
        "effect_size":   [  1,   1,  -1,   1,  -1,  -1,  -2,   1],
        "trait_id":      [  0,   0,   0,   0,   0,   0,   0,   0],
        "causal_allele": ["1", "1", "1", "1", "1", "1", "1", "1"]}
trait_df = pd.DataFrame(data)
trait_df
```

We obtain **edge effects** by using the {py:func}`edge_effect` function:

```{code-cell}
edge_effect_df = tstrait.edge_effect(ts, trait_df)
edge_effect_df
```

On edge 6, from node 7 to node 4, mutations 0 and 3 introduce causal alleles,
each with effect `1`, so the edge effect is `1 + 1 = 2`.
Edge 0 contains no mutations, so its edge effect is `0`.

TODO: Gregor to check these IDs

We obtain **edge genetic values** by using the `level="edge"` argument
of the {py:func}`genetic_value` function:

```{code-cell}
edge_value_df = tstrait.genetic_value(ts, trait_df, level="edge")
edge_value_df
```

Edge 1 inherits the value `2` from edge 6 over `[0, 5)` and
has its own effect of `-1`, so its edge genetic value is `2 + (-1) = 1`.

TODO: Gregor to check these IDs

We obtain **node genetic values** by using the `level="node"` argument
of the {py:func}`genetic_value` function:

```{code-cell}
node_value_df = tstrait.genetic_value(ts, trait_df, level="node")
node_value_df
```

Note that the node 1 is recombinant:
it inherits
edge 1 over `[0, 5)`, with genetic value `1` and
edge 2 over `[5, 10)`, with genetic value `-1`.
Its node genetic value is therefore `1 + (-1) = 0`.

We obtain **individual genetic values** by using the `level="individual"` argument
of the {py:func}`genetic_value` function, which is the default level:

```{code-cell}
individual_value_df = tstrait.genetic_value(ts, trait_df)
individual_value_df
```

The individual genetic values are `[2, -4]`.
Individual 0 contains nodes 0 and 1, so its value is `2 + 0 = 2`.
Individual 1 contains nodes 2 and 3, so its value is `-2 + (-2) = -4`.

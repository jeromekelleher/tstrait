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
this page will describe relationships between individual, node, and edge genetic values, and underlying allele effects.

**Learning Objectives**

After this page, you will be able to:

- Understand how to obtain genetic value in tstrait for individuals, nodes, and edges
- Understand how allele effects give rise to different genetic values
- Understand how to use allele effect model and mutation effect model in tstrait (TODO: later - move to a separate page?)

# Algorithm Overview

The tstrait algorithm used for computing
{ref}`individuals' genetic values<genetic_value_doc>` from allele effects
can also be used to compute node and edge genetic values.
These are related as follows:

- *edge effect* is a sum of allele effects on the edge.
- *edge "genetic" value* is a sum of edge effects above the edge and the edge.
- *node "genetic" value* is a sum of "genetic" values of incoming edges into the node.
- *individual genetic value* is a sum of "genetic" values of nodes of the individual.

Above we have quoted the term "genetic" when refering to nodes and edges,
because that term is usually used for individuals.
For consistency across different tstrait outputs,
we use the same term (genetic value) for individuals, nodes, and edges.

# Example

We will demonstrate the concepts of edge, node, and individual genetic values with
a tiny example so we can follow the computations.

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
# 0 [1, 0, 0, 0, 1, 0, 0, 0, 0, 0]
# 1 [1, 0, 1, 0, 1, 0, 1, 1, 0, 1]
# 2 [0, 0, 0, 0, 0, 0, 1, 1, 0, 0]
# 3 [0, 1, 0, 0, 0, 0, 1, 0, 1, 0]

#        mutation           0    1    2    3    4    5    6    7
data = {"site_id":       [  0,   1,   2,   4,   6,   7,   8,   9],
        "effect_size":   [  2,   1,  -1,  -1,  -1,  -1,   1,   1],
        "trait_id":      [  0,   0,   0,   0,   0,   0,   0,   0],
        "causal_allele": ["1", "1", "1", "1", "1", "1", "1", "1"]}
trait_df = pd.DataFrame(data)
trait_df
```

We obtain **edge effects** by using the {py:func}`edge_effect` function:

```{code-cell}
edge_effect_df = tstrait.edge_effect(ts, trait_df)
edge_effect_df
# ALL VALUES ARE CORRECT!
```

The above output is a {py:class}`pandas.DataFrame` object with the following columns:

> - **trait_id**: Trait ID that will be used in multi-trait simulation.
> - **edge_id**: Edge ID inside the tree sequence input.
> - **effect_size**: Simulated edge effect.

We obtain **edge genetic values** by using the `mode="edge"` argument
of the {py:func}`genetic_value` function:

```{code-cell}
edge_value_df = tstrait.genetic_value(ts, trait_df, mode="edge")
edge_value_df
# ERROR: Edge 1 should have value 1, but has 0!
```

The above output is a {py:class}`pandas.DataFrame` object with the following columns:

> - **trait_id**: Trait ID that will be used in multi-trait simulation.
> - **edge_id**: Edge ID inside the tree sequence input.
> - **genetic_value**: Simulated edge genetic value.

We obtain **node genetic values** by using the `mode="node"` argument
of the {py:func}`genetic_value` function:

```{code-cell}
node_value_df = tstrait.genetic_value(ts, trait_df, mode="node")
node_value_df
# ERROR: Node 1 should have value 0, but has -1 (due to incoming edge 1 having wrong value)
```

The above output is a {py:class}`pandas.DataFrame` object with the following columns:

> - **trait_id**: Trait ID that will be used in multi-trait simulation.
> - **node_id**: Node ID inside the tree sequence input.
> - **genetic_value**: Simulated node genetic value.

We obtain **individual genetic values** by using the `mode="individual"` argument
of the {py:func}`genetic_value` function, which is the default mode:

```{code-cell}
individual_value_df = tstrait.genetic_value(ts, trait_df)
individual_value_df
# ERROR: Individual 0 should have value 1, but has 0 (due to incoming edge 1 / node 1 having wrong value)!
```

The above output is a {py:class}`pandas.DataFrame` object with the following columns:

> - **trait_id**: Trait ID that will be used in multi-trait simulation.
> - **individual_id**: Individual ID inside the tree sequence input.
> - **genetic_value**: Simulated Individual genetic value.

# Allele and mutation effect models

TODO: Later
Show how this is done in tstrait.
Move to a separate page?

tstrait by default uses the allele effect model.
This model assumes that an allele has appeared with mutation once
in the tree sequence and that the effect of this allele is the same
for all individuals carrying this allele.

tstrait also supports the mutation effect model.
This models assumes that an allele can appear with mutation multiple times
in the tree sequence and that the effect of the resulting alleles differs
between mutation events and correspondingly for individuals carrying these alleles.
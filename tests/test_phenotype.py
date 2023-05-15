"""
import pytest

import msprime
import numpy as np
import tskit
import pandas as pd

import tstrait.phenotype as sim_pheno

class Test_update_node_values_array_access:
    def test_update(self):
        # Construct an artificial tree with 6 leaves with node 7 and 9 missing from the tree
        # We will be having mutation in every edge having genetic value 1 to be sure that all the edges are captured in the tree traversal algorithm
        root = 12
        left_child_array = np.array([-1,-1,-1,-1,-1,-1,0,-1,2,-1,4,6,11,12])
        right_sib_array = np.array([1,-1,3,-1,5,-1,8,-1,-1,-1,-1,10,-1,-1])
        node_values = np.array([1,1,1,1,1,1,1,0,1,0,1,1,1])
        
        # Create an arbitrary genotype vector
        rng = np.random.default_rng(1)
        G = rng.normal(loc=0, scale=1, size=6)        
        G1 = np.copy(G)
        sim_pheno.update_node_values_array_access(root, left_child_array, right_sib_array, node_values, G)
        
        actual_G = np.array([4,4,4,4,3,3]) + G1
        
        assert np.array_equal(G, actual_G)
      
    def binary_tree(self):
        # 3.00   6     
        #     ┊ ┏━┻━┓    ┊
        # 2.00┊ ┃   5   ┊
        #     ┊ ┃ ┏━┻┓   ┊
        # 1.00┊ ┃ ┃  4  ┊
        #     ┊ ┃ ┃ ┏┻┓  ┊
        # 0.00 0 1 2 3 
        #              
        # Mutation in all edges
        ts = tskit.Tree.generate_comb(4, span=10).tree_sequence
        tables = ts.dump_tables()
        for j in range(6):
            tables.sites.add_row(j, "A")
            tables.mutations.add_row(site=j, node=j, derived_state="T")
        ts = tables.tree_sequence()
        tree = ts.first()
        node_values = np.zeros(ts.num_nodes)
        for mut in tree.mutations():
            node_values[mut.node] += 1
        
        rng = np.random.default_rng(2)
        G = rng.normal(loc=0, scale=1, size=ts.num_samples)
        G1 = np.copy(G)
        sim_pheno.update_node_values_array_access(tree.root, tree.left_child_array, tree.right_sib_array, node_values, G)
        
        actual_G = np.array([1,2,3,3]) + G1
        
        assert np.array_equal(G, actual_G)
    
    def non_binary_tree(self):

        # 2.00      7   
        #     ┊ ┏━┏━━┏━┻━━━┓   ┊
        # 1.00┊ ┃ ┃ ┃    6   ┊
        #     ┊ ┃ ┃ ┃  ┏━┳┻┳━┓ ┊
        # 0.00 0 1 2  3 4 5
        #              
        # Mutation in all edges
        ts = tskit.Tree.generate_balanced(6, arity=4, span=10).tree_sequence
        tables = ts.dump_tables()
        for j in range(7):
            tables.sites.add_row(j, "A")
            tables.mutations.add_row(site=j, node=j, derived_state="T")
        ts = tables.tree_sequence()
        tree = ts.first()
        node_values = np.zeros(ts.num_nodes)
        for mut in tree.mutations():
            node_values[mut.node] += 1
        
        rng = np.random.default_rng(3)
        G = rng.normal(loc=0, scale=1, size=ts.num_samples)
        G1 = np.copy(G)
        sim_pheno.update_node_values_array_access(tree.root, tree.left_child_array, tree.right_sib_array, node_values, G)
        
        actual_G = np.array([1,1,1,2,2,2]) + G1
        
        assert np.array_equal(G, actual_G)

class Test_environment:
    @pytest.mark.parametrize("size", [1, 2, 10, 100])
    @pytest.mark.parametrize("h2", [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1])
    def test_environment(self, size, h2):
        rng = np.random.default_rng(1)
        G = rng.normal(size = size)
        phenotype, E = sim_pheno.environment(G, h2, 1, rng)
        assert len(phenotype) == size
        assert len(E) == size
        assert np.issubdtype(phenotype.dtype, np.floating)
        assert np.issubdtype(E.dtype, np.floating)
    
    @pytest.mark.parametrize("h2", [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1])
    def test_environment_G_error(self, h2):
        rng = np.random.default_rng(2)
        G = rng.normal(size = 0)
        with pytest.raises(ValueError, match="No individuals in the simulation model"):
            phenotype, E = sim_pheno.environment(G, h2, 1, rng)
    
    @pytest.mark.parametrize("size", [1, 2, 10, 100])    
    @pytest.mark.parametrize("h2", [-0.2, -0.01, 1.01, 1.2, "a", "1", [0.4,0.3]])
    def test_h2error(self, size, h2):
        G = np.random.normal(size = size)
        rng = np.random.default_rng(3)        
        with pytest.raises(ValueError, match="Heritability should be 0 <= h2 <= 1"):
            phenotype, E = sim_pheno.environment(G, h2, 1, rng)   

class Test_choose_causal:
    @pytest.mark.parametrize("num_mutations", [10, 100, 1000])
    @pytest.mark.parametrize("num_causal", [1, 2, 10])
    @pytest.mark.parametrize("trait_sd", [0.1, 0.3, 0.5, 0.7, 0.9])
    @pytest.mark.parametrize("trait_mean", [-1,0,1])
    def test_choose_causal(self, num_mutations, num_causal, trait_sd, trait_mean):
        rng = np.random.default_rng(1)
        mutation_id, beta = sim_pheno.choose_causal(num_mutations, num_causal, trait_mean, trait_sd, rng)
        assert len(beta) == num_causal
        assert len(mutation_id) == num_causal
        assert min(mutation_id) >= 0
        assert max(mutation_id) < num_mutations
        assert np.issubdtype(beta.dtype, np.floating)
        assert np.issubdtype(mutation_id.dtype, np.integer)

    @pytest.mark.parametrize("num_mutations", [10, 100, 1000])
    @pytest.mark.parametrize("trait_sd", [0.1, 0.3, 0.5, 0.7, 0.9])   
    @pytest.mark.parametrize("trait_mean", [-1,0,1])    
    def test_zero(self, num_mutations, trait_sd, trait_mean):
        rng = np.random.default_rng(1)
        num_causal = 0
        mutation_id, beta = sim_pheno.choose_causal(num_mutations, num_causal, trait_mean, trait_sd, rng)
        assert len(beta) == 0
        assert len(mutation_id) == 0
    
    @pytest.mark.parametrize("addition", [1, 2, 10, 100])
    @pytest.mark.parametrize("num_mutations", [1, 2, 10, 100])
    @pytest.mark.parametrize("trait_sd", [0.1, 0.3, 0.5, 0.7, 0.9]) 
    @pytest.mark.parametrize("trait_mean", [-1,0,1])
    def test_error_num_sites(self, addition, num_mutations, trait_sd, trait_mean):
        rng = np.random.default_rng(2)
        with pytest.raises(ValueError, match="There are more causal sites than the number of mutations inside the tree sequence"):
            num_causal = num_mutations + addition
            mutation_id, beta = sim_pheno.choose_causal(num_mutations, num_causal, trait_mean, trait_sd, rng)
    
    @pytest.mark.parametrize("num_mutations", [10, 100, 1000])
    @pytest.mark.parametrize("trait_sd", [0.1, 0.3, 0.5, 0.7, 0.9]) 
    @pytest.mark.parametrize("num_causal", [-1, 0]) 
    @pytest.mark.parametrize("trait_mean", [-1,0,1])
    def test_error_num_causal(self, num_mutations, trait_sd, num_causal, trait_mean):
        rng = np.random.default_rng(3)
        with pytest.raises(ValueError, match = "Number of causal sites should be a non-negative integer"):
            mutation_id, beta = sim_pheno.choose_causal(num_mutations, num_causal, trait_mean, trait_sd, rng)
        
    @pytest.mark.parametrize("num_mutations", [10, 100, 1000])
    @pytest.mark.parametrize("trait_sd", [0, -0.1, -1, -10])
    @pytest.mark.parametrize("num_causal", [1, 2, 10]) 
    @pytest.mark.parametrize("trait_mean", [-1,0,1])    
    def test_error_trait_sd(self, num_mutations, trait_sd, num_causal, trait_mean):
        rng = np.random.default_rng(np.random.randint(5))
        with pytest.raises(ValueError, match = "Standard deviation should be a non-negative number"):
            mutation_id, beta = sim_pheno.choose_causal(num_mutations, num_causal, trait_mean, trait_sd, rng)



class Test_genetic_value:
    @pytest.mark.parametrize("seed", [1, 3, 5, 7, 9])
    @pytest.mark.parametrize("num_causal", [1, 2, 10])
    def test_genetic_value(self, seed, num_causal):
        rng = np.random.default_rng(seed)
        ts = msprime.sim_ancestry(100, sequence_length=10_000, recombination_rate=1e-8,population_size=10**4, random_seed=seed)
        ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed, model="binary", discrete_genome=False)
        mutation_id, beta = sim_pheno.choose_causal(ts.num_mutations, num_causal,0, 1, rng)
        G, location, mutation_list = sim_pheno.genetic_value(ts, mutation_id, beta)
        assert len(G) == ts.num_individuals
        assert len(location) == num_causal
        assert len(mutation_list) == num_causal
        assert np.issubdtype(G.dtype, np.floating)
        assert np.issubdtype(location.dtype, np.integer) or np.issubdtype(location.dtype, np.floating)
        assert np.issubdtype(mutation_list.dtype, np.integer)
    
    def test_tree_sequence_binary(self):
        ts = tskit.Tree.generate_comb(6, span=20).tree_sequence
        tables = ts.dump_tables()
        for j in range(10):
            tables.sites.add_row(j, "A")
            tables.mutations.add_row(site=j, node=j, derived_state="T")
        ts = tables.tree_sequence()
        mutation_id = np.array(list(range(10)))
        beta = np.array(list(range(10)))
        G = sim_pheno.genetic_value(ts, mutation_id, beta)
        
        G = G[::2] + G[1::2]
        
        G_actual = np.array([10, 46, 69])

        assert np.array_equal(G, G_actual)
        
    def test_tree_sequence_nonbinary(self):
        rng = np.random.default_rng(2)
        ts = tskit.Tree.generate_balanced(6, arity=4, span=10).tree_sequence
        tables = ts.dump_tables()
        for j in range(7):
            tables.sites.add_row(j, "A")
            tables.mutations.add_row(site=j, node=j, derived_state="T")
        ts = tables.tree_sequence()
        mutation_id = np.array([0, 1, 3, 6])
        beta = np.repeat(1,4)
        G = sim_pheno.genetic_value(ts, mutation_id, beta)
        G = G[::2] + G[1::2]
        G_actual = np.array([2,2,2])
        
        assert np.array_equal(G, G_actual)
              
class Test_phenotype_sim:
    @pytest.mark.parametrize("seed", [1, 2, 3, 4])
    @pytest.mark.parametrize("num_causal", [1, 2, 10])
    def test_phenotype_sim(self, seed, num_causal):
        size = 20
        ts = msprime.sim_ancestry(size, sequence_length=100_000, recombination_rate=1e-8,population_size=10**4, random_seed=seed)
        ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed, model="binary", discrete_genome=False)
        pheno_df, gene_df = sim_pheno.phenotype_sim(ts, num_causal, trait_sd=1, h2=0.3, seed=seed)
        assert pheno_df.isnull().values.any() == False
        assert gene_df.isnull().values.any() == False
        assert len(pheno_df) == ts.num_individuals
        assert len(gene_df) == num_causal
    
    @pytest.mark.parametrize("seed", [1, 2, 3, 4])
    @pytest.mark.parametrize("num_causal", [1, 2, 10])
    def test_no_mutation(self, seed, num_causal):
        size = 20
        rng = np.random.default_rng(seed)
        ts = msprime.sim_ancestry(size, sequence_length=100_000, recombination_rate=1e-8,population_size=10**4, random_seed=seed)
        with pytest.raises(ValueError, match = "No mutation in the provided data"):
            pheno_df, gene_df = sim_pheno.phenotype_sim(ts, num_causal, trait_sd=1, h2=0.3, seed=seed)
"""
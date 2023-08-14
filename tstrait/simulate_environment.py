import numbers

import numpy as np

from .base import _check_dataframe


class EnvSimulator:
    """Simulator class to simulate environmental noise of individuals.

    :param genetic_df: Pandas dataframe that includes genetic value of individuals.
        It must include individual ID, genetic value and trait ID.
    :type genetic_df: pandas.DataFrame
    :param h2: Narrow-sense heritability, which will be used to simulate environmental
        noise. Narrow-sense heritability must be between 0 and 1.
    :type h2: float or list or numpy.ndarray
    :param random_seed: The random seed. If this is not specified or None, simulation
        will be done randomly.
    :type random_seed: int
    """

    def __init__(self, genetic_df, h2, random_seed):
        self.genetic_df = genetic_df[["trait_id", "individual_id", "genetic_value"]]
        self.h2 = h2
        self.rng = np.random.default_rng(random_seed)

    def _sim_env(self, var, h2):
        """Simulate environmental noise based on variance and narrow-sense
        heritability
        """
        env_std = np.sqrt((1 - h2) / h2 * var)
        env_noise = self.rng.normal(loc=0.0, scale=env_std)

        return env_noise

    def sim_environment(self):
        """Simulate environmental values based on genetic values of individuals and
        narrow-sense heritability
        """
        df = self.genetic_df.copy()
        h2_array = np.take(self.h2, self.genetic_df.trait_id)

        grouped = df.groupby("trait_id", sort=False)["genetic_value"]
        var_array = grouped.transform("var")

        df["environmental_noise"] = self._sim_env(var_array, h2_array)

        df["phenotype"] = df["genetic_value"] + df["environmental_noise"]

        return df


def sim_env(genetic_df, h2, random_seed=None):
    """Simulates environmental noise and obtain phenotype.

    :param genetic_df: Pandas dataframe that includes genetic value of individuals.
        It must include individual ID, genetic value and trait ID.
    :type genetic_df: pandas.DataFrame
    :param h2: Narrow-sense heritability, which will be used to simulate environmental
        noise. Narrow-sense heritability must be between 0 and 1.
    :type h2: float or list or numpy.ndarray
    :param random_seed: The random seed. If this is not specified or None, simulation
        will be done randomly.
    :type random_seed: int
    :return: Returns a pandas dataframe that includes individual ID, genetic value,
        environmental noise, phenotype, and trait ID.
    :rtype: pandas.DataFrame
    """
    genetic_df = _check_dataframe(
        genetic_df, ["trait_id", "individual_id", "genetic_value"], "genetic_df"
    )

    trait_id = genetic_df["trait_id"].unique()

    if np.min(trait_id) != 0 or np.max(trait_id) != len(trait_id) - 1:
        raise ValueError("trait_id must be consecutive and start from 0")

    if isinstance(h2, numbers.Real):
        h2 = np.array([h2])

    if len(h2) != len(trait_id):
        raise ValueError("Length of h2 must match the number of traits")

    if np.min(h2) <= 0 or np.max(h2) > 1:
        raise ValueError("Narrow-sense heritability must be 0 < h2 <= 1")

    simulator = EnvSimulator(
        genetic_df=genetic_df,
        h2=h2,
        random_seed=random_seed,
    )

    phenotype_df = simulator.sim_environment()

    return phenotype_df
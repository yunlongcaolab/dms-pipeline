# from J. Bloom Lab. https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html
from dms_variants.globalepistasis import (
    MonotonicSplineEpistasisGaussianLikelihood,
    MonotonicSplineEpistasisCauchyLikelihood,
    NoEpistasisGaussianLikelihood,
)
import jax.numpy as jnp
import numpy as np
import pandas as pd
import jax
from jax import grad, jit


def get_mutations(variants: list[str], split: str = " "):
    r"""
    Get the list of mutations from the given list of variants.
    """
    mutations = set()
    for variant in variants:
        for mut in variant.split(split):
            mutations.add(mut)
    return list(mutations)


@jit
def _calc_forward(data, params):
    r"""
    Compute the forward pass of the model.
    """
    linear_res = jnp.dot(data, params[0:-2]) + params[-2]
    return (jax.nn.sigmoid(linear_res) - jax.nn.sigmoid(params[-2])) * jnp.exp(
        params[-1]
    )


@jit
def _calc_forward_exp(data, params):
    r"""
    Compute the forward pass of the model.
    """
    linear_res = jnp.dot(data, jnp.exp(params[0:-2])) + params[-2]
    return (jax.nn.sigmoid(linear_res) - jax.nn.sigmoid(params[-2])) * jnp.exp(
        params[-1]
    )


@jit
def mse_loss_with_lasso(y_true: jnp.ndarray, y_pred: jnp.ndarray, lasso: float):
    return jnp.mean(jnp.square(y_true - y_pred)) + lasso * jnp.mean(jnp.abs(y_pred))
    # return jnp.mean(jnp.square(y_true - y_pred)) + lasso * jnp.mean(jnp.abs(params[:-2]))


# @jit
# def mse_loss_with_lasso_exp(y_true: jnp.ndarray, y_pred: jnp.ndarray, params: jnp.ndarray, lasso:float):
#     return jnp.mean(jnp.square(y_true - y_pred)) + lasso * jnp.mean(jnp.abs(y_pred))
#     # return jnp.mean(jnp.square(y_true - y_pred)) + lasso * jnp.mean(jnp.exp(params[:-2]))


class SigmoidEpistasis:
    def __init__(
        self,
        mutations: list[str],
        seed: int = 42,
        exp: bool = False,
        lasso: float = 0.0,
    ):
        self.mutations = set(mutations)
        self.randkey = jax.random.PRNGKey(seed)

        self.mutations = list(self.mutations)
        self.mut2idx = {mut: i for i, mut in enumerate(self.mutations)}

        self.exp = exp

        # initialize a parameter for each mutation (normal)
        self.params = jax.random.normal(self.randkey, (len(self.mutations),))
        self.final_mapping = jax.random.normal(
            self.randkey, (2,)
        )  # slope and intercept

        self.lasso = lasso

    def latent_to_observed(self, value):
        observed = jax.nn.sigmoid(value + self.final_mapping[0]) - jax.nn.sigmoid(
            self.final_mapping[0]
        )
        return observed * jnp.exp(self.final_mapping[1])

    def get_effects(
        self, cut_min: float = None, cut_max: float = None, normalize: bool = False
    ):
        r"""
        Get the effects of the mutations.
        """

        _params = self.params if not self.exp else jnp.exp(self.params)
        # effects = [self.latent_to_observed(_params[i]) for i in range(len(self.mutations))]
        effects = self.latent_to_observed(_params)

        if normalize:
            effects = (effects - jnp.min(effects)) / jnp.max(effects)

        else:
            if cut_min is not None:
                effects = jnp.maximum(effects, cut_min)
            if cut_max is not None:
                effects = jnp.minimum(effects, cut_max)

        return pd.DataFrame(
            {
                "mutation": self.mutations,
                "latent": _params.tolist(),
                "observed": effects,
            }
        )

    def construct_data(self, variants: dict[str, float], split=" "):
        r"""
        Construct the data matrix for the model.
        variants: A dictionary where keys are the variants and values are the corresponding values.
        A vector with value {0,1} where 1 indicates the mutation is present in the variant.
        """
        data = np.zeros((len(variants), len(self.mutations)))
        labels = np.zeros((len(variants),))
        for i, key in enumerate(variants.keys()):
            val = variants[key]
            for mut in key.split(split):
                data[i, self.mut2idx[mut]] = 1
            labels[i] = val

        return jnp.array(data), jnp.array(labels)

    def __call__(self, data: jnp.ndarray):
        r"""
        Compute the predicted values for the given data.
        data: [n_samples, n_total_mutations]
        """

        return (
            _calc_forward(data, self.get_all_params())
            if not self.exp
            else _calc_forward_exp(data, self.get_all_params())
        )

    def get_all_params(self):
        return jnp.concatenate([self.params, self.final_mapping])

    def set_all_params(self, params: jnp.ndarray):
        self.params = params[:-2]
        self.final_mapping = params[-2:]

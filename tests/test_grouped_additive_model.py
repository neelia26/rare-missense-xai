import unittest

import numpy as np
import torch

from src.modeling.grouped_additive import (
    RECONSTRUCTION_TOLERANCE,
    GroupedAdditiveNetwork,
    GroupedAdditiveSettings,
    _predict,
    _train_epochs,
    group_indices,
    set_deterministic_seed,
)


class GroupedAdditiveArchitectureTests(unittest.TestCase):
    def setUp(self):
        self.indexes = {"conservation": (0, 1), "constraint": (2,), "amino": (3, 4), "protein": (5,)}
        self.settings = GroupedAdditiveSettings(hidden_width=4)

    def test_group_outputs_reconstruct_logit_exactly(self):
        set_deterministic_seed(7)
        model = GroupedAdditiveNetwork(self.indexes, self.settings)
        values = torch.tensor(np.arange(18).reshape(3, 6), dtype=torch.float32)
        logits, contributions = model(values)
        reconstructed = model.intercept + contributions.sum(dim=1)
        torch.testing.assert_close(logits, reconstructed, rtol=0, atol=1e-7)
        self.assertLess(RECONSTRUCTION_TOLERANCE, 1e-4)

    def test_each_group_has_one_scalar_output(self):
        model = GroupedAdditiveNetwork(self.indexes, self.settings)
        _, contributions = model(torch.zeros((5, 6)))
        self.assertEqual(tuple(contributions.shape), (5, 4))

    def test_group_index_resolution_is_complete_and_disjoint(self):
        names = ["a", "b", "c"]
        resolved = group_indices(names, {"one": ("a", "c"), "two": ("b",)})
        self.assertEqual(resolved, {"one": (0, 2), "two": (1,)})
        with self.assertRaises(ValueError):
            group_indices(names, {"one": ("a", "b"), "two": ("b", "c")})

    def test_seed_reproduces_initial_parameters(self):
        set_deterministic_seed(19)
        first = GroupedAdditiveNetwork(self.indexes, self.settings)
        first_state = {key: value.detach().clone() for key, value in first.state_dict().items()}
        set_deterministic_seed(19)
        second = GroupedAdditiveNetwork(self.indexes, self.settings)
        for key, value in second.state_dict().items():
            torch.testing.assert_close(first_state[key], value, rtol=0, atol=0)

    def test_same_seed_reproduces_tiny_training_predictions(self):
        values = np.arange(72, dtype=np.float32).reshape(12, 6) / 20
        labels = np.array([0, 1] * 6, dtype=np.float32)
        predictions = []
        for _ in range(2):
            set_deterministic_seed(31)
            model = GroupedAdditiveNetwork(self.indexes, self.settings)
            _train_epochs(model, values, labels, self.settings, seed=31, epochs=3)
            predictions.append(_predict(model, values)[1])
        np.testing.assert_array_equal(predictions[0], predictions[1])


if __name__ == "__main__":
    unittest.main()

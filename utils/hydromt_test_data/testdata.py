"""Custom HydroMT steps for the Wflow test data, as a model component.

The HydroMT workflow files in ``hydromt/`` only contain HydroMT-Wflow steps and the steps
of this component, prefixed with ``testdata.``. The component is added to the model with
``model.add_component("testdata", TestDataComponent(model))``, so ``Model.build`` and
``Model.update`` resolve and validate these steps like any other step. The component has
no data of its own: the steps add maps, forcing and config entries to the HydroMT-Wflow
components. The implementation of the steps is in ``steps.py``.
"""

import functools
import inspect

from hydromt.model.components import ModelComponent
from hydromt.model.steps import hydromt_step

import steps


def _component_step(func):
    """Expose ``func(model, ...)`` as a HydroMT step ``method(self, ...)`` of the component."""

    @functools.wraps(func)
    def method(self, *args, **kwargs):
        return func(self.model, *args, **kwargs)

    parameters = list(inspect.signature(func).parameters.values())[1:]
    self_parameter = inspect.Parameter("self", inspect.Parameter.POSITIONAL_OR_KEYWORD)
    method.__signature__ = inspect.Signature([self_parameter, *parameters])
    return hydromt_step(method)


class TestDataComponent(ModelComponent):
    """Steps that adapt a HydroMT-Wflow model to the Wflow.jl test data."""

    __test__ = False  # not a pytest test class

    def read(self) -> None:
        """Nothing to read: the steps only modify other components."""

    def write(self) -> None:
        """Nothing to write: the steps only modify other components."""

    @hydromt_step
    def read_default_config(self) -> None:
        """Replace the config by the HydroMT-Wflow default config of the model type.

        Used when a sediment model is built as an update of an sbm model.
        """
        self.model.config.read(self.model.config._default_template_filename)

    setup_lst_parameters = _component_step(steps.lst_parameters)
    setup_lst_forcing = _component_step(steps.era5_lst_forcing)
    setup_groundwater_maps = _component_step(steps.groundwater_maps)
    setup_layered_map = _component_step(steps.stack_layers)
    setup_point_map = _component_step(steps.point_map)
    setup_point_timeseries = _component_step(steps.point_timeseries)

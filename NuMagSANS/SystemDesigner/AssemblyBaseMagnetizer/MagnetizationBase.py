"""Parameter bases for local magnetization-field materialization.

This module mirrors the parameter-distribution logic of AssemblyBaseDesigner,
but the semantics are magnetic-field templates instead of crystal templates.
It does not place objects in global space and it does not rotate whole objects.
"""

import numpy as np

from ..AssemblyBaseDesigner.AssemblyBase import (
    DISTRIBUTION_TYPES,
)


def _as_template_list(field_templates):
    """Normalize one field-template name or many names to a list."""
    if isinstance(field_templates, str):
        return [field_templates]

    templates = list(field_templates)
    if not templates:
        raise ValueError("At least one magnetization field template must be provided.")

    if not all(isinstance(template, str) for template in templates):
        raise TypeError("Magnetization field template names must be strings.")

    return templates


def _normalize_template_params(field_templates, field_template_params):
    """Normalize parameter declarations to template -> parameter list."""
    if len(field_templates) == 1 and isinstance(field_template_params, (list, tuple, set)):
        return {field_templates[0]: list(field_template_params)}

    if not isinstance(field_template_params, dict):
        raise TypeError(
            "field_template_params must be a parameter list for one template or "
            "a dictionary mapping field-template names to parameter lists."
        )

    if len(field_templates) == 1 and field_templates[0] not in field_template_params:
        return {field_templates[0]: list(field_template_params.keys())}

    normalized = {}
    for template in field_templates:
        if template not in field_template_params:
            raise KeyError(f"Missing parameter declaration for field template: {template}")

        params = list(field_template_params[template])
        if not params:
            raise ValueError(f"Field template {template} must define at least one parameter.")

        if not all(isinstance(param, str) for param in params):
            raise TypeError(f"Parameter names for field template {template} must be strings.")

        normalized[template] = params

    return normalized


def _validate_distribution_spec(spec, label):
    """Validate one distribution specification."""
    if "distribution" not in spec:
        raise KeyError(f"{label} is missing a distribution type.")

    distribution = spec["distribution"]
    if distribution not in DISTRIBUTION_TYPES:
        raise ValueError(f"{label} has unknown distribution {distribution!r}. Allowed values are {DISTRIBUTION_TYPES}.")

    required_keys = {
        "constant": ("value",),
        "uniform": ("low", "high"),
        "normal": ("mean", "sigma"),
        "lognormal": ("mean", "sigma", "mean_type"),
    }

    missing_keys = [key for key in required_keys[distribution] if key not in spec]
    if missing_keys:
        raise KeyError(f"{label} is missing distribution keys: {missing_keys}")


def _normalize_distribution_props(field_templates, template_params, param_dist_props):
    """Normalize distribution declarations to template -> parameter -> spec."""
    if not isinstance(param_dist_props, dict):
        raise TypeError("param_dist_props must be a dictionary.")

    if len(field_templates) == 1 and not any(template in param_dist_props for template in field_templates):
        param_dist_props = {field_templates[0]: param_dist_props}

    normalized = {}
    for template in field_templates:
        if template not in param_dist_props:
            raise KeyError(f"Missing distribution properties for field template: {template}")

        normalized[template] = {}
        for param in template_params[template]:
            if param not in param_dist_props[template]:
                raise KeyError(f"Missing distribution properties for parameter {template}.{param}")

            spec = dict(param_dist_props[template][param])
            _validate_distribution_spec(spec, f"{template}.{param}")
            normalized[template][param] = spec

    return normalized


def magnetization_base(field_templates, field_template_params, param_dist_props, name=None):
    """Create a parameter-level base for local magnetization generation.

    The magnetization base defines which local vector-field templates can be
    used and how their field parameters are sampled. It does not define global
    object positions or global object rotations.
    """
    field_templates = _as_template_list(field_templates)
    template_params = _normalize_template_params(field_templates, field_template_params)
    distribution_props = _normalize_distribution_props(
        field_templates,
        template_params,
        param_dist_props,
    )

    return {
        "type": "magnetization_base",
        "name": name,
        "field_templates": {
            template: {
                "parameters": {param: distribution_props[template][param] for param in template_params[template]}
            }
            for template in field_templates
        },
    }


def _sample_distribution(spec, rng, size=None):
    """Sample values from one distribution specification."""
    distribution = spec["distribution"]

    if distribution == "constant":
        value = spec["value"]
        if size is None:
            return value
        return np.full(size, value)

    if distribution == "uniform":
        return rng.uniform(spec["low"], spec["high"], size=size)

    if distribution == "normal":
        return rng.normal(spec["mean"], spec["sigma"], size=size)

    if distribution == "lognormal":
        mean = spec["mean"]
        sigma = spec["sigma"]
        if spec["mean_type"] == "arithmetic":
            mean = np.log(mean) - 0.5 * sigma**2
        return rng.lognormal(mean=mean, sigma=sigma, size=size)

    raise ValueError(f"Unknown distribution type: {distribution}")


def sample_field_params(magnetization_base_data, field_template_name, size=None, seed=None):
    """Sample concrete vector-field parameters from a magnetization base."""
    if magnetization_base_data.get("type") != "magnetization_base":
        raise ValueError("magnetization_base_data must be created by magnetization_base().")

    templates = magnetization_base_data["field_templates"]
    if field_template_name not in templates:
        raise KeyError(f"Unknown field template: {field_template_name}")

    rng = np.random.default_rng(seed)
    return {
        param: _sample_distribution(spec, rng, size=size)
        for param, spec in templates[field_template_name]["parameters"].items()
    }

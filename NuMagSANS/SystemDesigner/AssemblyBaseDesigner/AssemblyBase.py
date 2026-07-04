import numpy as np

DISTRIBUTION_TYPES = ("constant", "uniform", "normal", "lognormal")


def constant(value):
    """Define a template parameter that always has the same value."""
    return {
        "distribution": "constant",
        "value": value,
    }


def uniform(low, high):
    """Define a template parameter sampled uniformly from [low, high)."""
    if high <= low:
        raise ValueError("uniform distribution requires high > low.")

    return {
        "distribution": "uniform",
        "low": low,
        "high": high,
    }


def normal(mean, sigma):
    """Define a template parameter sampled from a normal distribution."""
    if sigma < 0:
        raise ValueError("normal distribution requires sigma >= 0.")

    return {
        "distribution": "normal",
        "mean": mean,
        "sigma": sigma,
    }


def lognormal(mean, sigma, mean_type="arithmetic"):
    """Define a template parameter sampled from a log-normal distribution.

    Parameters
    ----------
    mean : float
        Mean radius or log-space mean, depending on ``mean_type``.
    sigma : float
        Standard deviation in log space.
    mean_type : {"arithmetic", "log"}, optional
        If ``"arithmetic"`` (default), ``mean`` is interpreted as the desired
        arithmetic mean of the sampled values. If ``"log"``, ``mean`` is used
        directly as the NumPy log-normal ``mean`` parameter.
    """
    if mean <= 0 and mean_type == "arithmetic":
        raise ValueError("arithmetic lognormal mean must be > 0.")

    if sigma < 0:
        raise ValueError("lognormal distribution requires sigma >= 0.")

    if mean_type not in ("arithmetic", "log"):
        raise ValueError('mean_type must be either "arithmetic" or "log".')

    return {
        "distribution": "lognormal",
        "mean": mean,
        "sigma": sigma,
        "mean_type": mean_type,
    }


def _as_template_list(ct_temps):
    """Normalize one template name or many template names to a list."""
    if isinstance(ct_temps, str):
        return [ct_temps]

    templates = list(ct_temps)
    if not templates:
        raise ValueError("At least one crystal template must be provided.")

    if not all(isinstance(template, str) for template in templates):
        raise TypeError("Crystal template names must be strings.")

    return templates


def _normalize_template_params(ct_temps, ct_temps_params):
    """Normalize template parameter declarations to template -> parameter list."""
    if len(ct_temps) == 1 and isinstance(ct_temps_params, (list, tuple, set)):
        return {ct_temps[0]: list(ct_temps_params)}

    if not isinstance(ct_temps_params, dict):
        raise TypeError(
            "ct_temps_params must be a parameter list for one template or a "
            "dictionary mapping template names to parameter lists."
        )

    if len(ct_temps) == 1 and ct_temps[0] not in ct_temps_params:
        return {ct_temps[0]: list(ct_temps_params.keys())}

    normalized = {}
    for template in ct_temps:
        if template not in ct_temps_params:
            raise KeyError(f"Missing parameter declaration for template: {template}")

        params = list(ct_temps_params[template])
        if not params:
            raise ValueError(f"Template {template} must define at least one parameter.")

        if not all(isinstance(param, str) for param in params):
            raise TypeError(f"Parameter names for template {template} must be strings.")

        normalized[template] = params

    return normalized


def _normalize_distribution_props(ct_temps, template_params, param_dist_props):
    """Normalize parameter distribution declarations to template -> parameter -> spec."""
    if not isinstance(param_dist_props, dict):
        raise TypeError("param_dist_props must be a dictionary.")

    if len(ct_temps) == 1 and not any(template in param_dist_props for template in ct_temps):
        param_dist_props = {ct_temps[0]: param_dist_props}

    normalized = {}
    for template in ct_temps:
        if template not in param_dist_props:
            raise KeyError(f"Missing distribution properties for template: {template}")

        normalized[template] = {}
        for param in template_params[template]:
            if param not in param_dist_props[template]:
                raise KeyError(f"Missing distribution properties for parameter {template}.{param}")

            spec = dict(param_dist_props[template][param])
            _validate_distribution_spec(spec, f"{template}.{param}")
            normalized[template][param] = spec

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


def assembly_base(ct_temps, ct_temps_params, param_dist_props, name=None):
    """Create a parameter-level base for later assembly generation.

    The assembly base does not place objects in space. It only defines which
    crystal templates can be used and how their template parameters are sampled.
    Translation and rotation are added later by an assembly organizer.

    Parameters
    ----------
    ct_temps : str or iterable of str
        Crystal template name or names, for example ``"sc_sphere_crystal"``.
    ct_temps_params : iterable or dict
        For one template, this can be a list of parameter names. For multiple
        templates, pass a dictionary mapping template names to parameter lists.
    param_dist_props : dict
        Distribution specification for each template parameter. For one template
        this can directly map parameter names to distribution specs. For
        multiple templates, pass a nested dictionary.
    name : str, optional
        Human-readable name for this assembly base.

    Returns
    -------
    dict
        Assembly base dictionary containing template names and validated
        parameter distribution specifications.
    """
    ct_temps = _as_template_list(ct_temps)
    template_params = _normalize_template_params(ct_temps, ct_temps_params)
    distribution_props = _normalize_distribution_props(
        ct_temps,
        template_params,
        param_dist_props,
    )

    return {
        "type": "assembly_base",
        "name": name,
        "crystal_templates": {
            template: {
                "parameters": {param: distribution_props[template][param] for param in template_params[template]}
            }
            for template in ct_temps
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


def sample_template_params(assembly_base_data, template_name, size=None, seed=None):
    """Sample concrete crystal-template parameters from an assembly base.

    Parameters
    ----------
    assembly_base_data : dict
        Dictionary returned by ``assembly_base``.
    template_name : str
        Name of the crystal template whose parameters should be sampled.
    size : int, optional
        Number of parameter sets to sample. If omitted, scalar values are
        returned. If provided, each parameter value is an array of length
        ``size``.
    seed : int, optional
        Random seed for reproducible sampling.

    Returns
    -------
    dict
        Sampled template parameter values.
    """
    if assembly_base_data.get("type") != "assembly_base":
        raise ValueError("assembly_base_data must be created by assembly_base().")

    templates = assembly_base_data["crystal_templates"]
    if template_name not in templates:
        raise KeyError(f"Unknown crystal template: {template_name}")

    rng = np.random.default_rng(seed)
    return {
        param: _sample_distribution(spec, rng, size=size)
        for param, spec in templates[template_name]["parameters"].items()
    }

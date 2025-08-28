"""
inputs.py

A tiny, dependency-light pattern to make inputs obvious, validated,
printable, and reproducible. No typing libs, no globals.

Usage in any script (e.g., calibrate.py):

    from calibrate_inputs_registry import (
        collect_inputs, validate_inputs, log_inputs_summary,
        inputs_markdown, save_inputs, InputsRegistry, compose_specs, pick,
    )

    SPEC_CALIBRATE = [
        {"path":"cli.standard_name","source":"args","key":"standard_name","required":True,"units":"name","desc":"Spectrophotometric standard"},
        {"path":"cli.output_folder","source":"args","key":"output_folder","required":True,"units":"path","desc":"Output dir"},
        {"path":"cli.extraction_radius","source":"args","key":"extraction_radius","required":False,"default":2.0,"units":"pix","desc":"Aperture radius","validate": lambda v: (v is None) or (_is_number(v) and v>0)},
        {"path":"dataset.reduced_files","source":"dataset","key":"reduced_files","required":True,"units":"list[str]","desc":"Reduced frames"},
        {"path":"config.fiber_radius","source":"config","key":"fiber_radius","required":True,"units":"pix","desc":"Fiber footprint radius","validate": lambda v: _is_number(v) and v>0},
    ]

    reg = InputsRegistry(SPEC_CALIBRATE)
    inputs = reg.collect(dataset_manifest, args, config_dict)
    reg.validate(inputs)
    reg.log(inputs)
    reg.save(inputs, Path(inputs["cli"]["output_folder"]) / "inputs.json")

Then pass only what is needed to downstream functions, e.g.:

    r = inputs["cli"]["extraction_radius"]
    files = inputs["dataset"]["reduced_files"]

"""

import json
from pathlib import Path
import logging

# Module logger (no external logger dependency)
logger = logging.getLogger('antigen.inputs')
# ---------- 1) Public, spec-agnostic API ----------

def collect_inputs(dataset_manifest, args, config_dict, spec):
    """Collect inputs from dataset/args/config into a nested dict per `spec`.

    Args:
        dataset_manifest (dict): Dictionary for dataset values.
        args: Namespace-like object (attributes for CLI values).
        config_dict (dict): Instrument configuration dictionary.
        spec (list[dict]): Field descriptors with keys like path/source/key/...

    Returns:
        dict: Nested mapping with top-level keys such as 'cli', 'dataset', 'config'.
    """
    inputs = {"cli": {}, "dataset": {}, "config": {}}

    for f in spec:
        path = f["path"]
        source = f["source"]
        key = f["key"]
        default = f.get("default")

        if source == "args":
            value = _safe_getattr(args, key, default)
        elif source == "dataset":
            value = (dataset_manifest or {}).get(key, default)
        elif source == "config":
            value = (config_dict or {}).get(key, default)
        else:
            value = default

        _nested_set(inputs, path, value)

    return inputs


def validate_inputs(inputs, spec):
    """Validate presence and sanity of inputs per `spec`.

    Raises:
        ValueError: with a readable list of problems.
    """
    problems = []

    for f in spec:
        path = f["path"]
        required = bool(f.get("required"))
        validate = f.get("validate")
        desc = f.get("desc", path)

        value = _nested_get(inputs, path)

        if required and _is_missing(value):
            problems.append("Missing required input: %s — %s" % (path, desc))
            continue

        if validate and not _call_validator(validate, value):
            problems.append("Invalid value for %s: %r" % (path, value))

    if problems:
        msg = "Input validation failed:\n" + "\n".join(problems)
        raise ValueError(msg)


def inputs_markdown(inputs, spec):
    """Render a markdown table documenting current inputs for a `spec`."""
    lines = [
        "| Key | Value | Units | Source | Required | Description |",
        "| --- | ----- | ----- | ------ | -------- | ----------- |",
    ]

    for f in spec:
        path = f["path"]
        units = f.get("units", "") or ""
        src = f["source"]
        req = "yes" if f.get("required") else "no"
        desc = f.get("desc", "")
        raw = _nested_get(inputs, path)
        value = _fmt_value(raw)
        lines.append("| %s | %s | %s | %s | %s | %s |" % (path, value, units, src, req, desc))

    return "\n".join(lines)


def log_inputs_summary(inputs, spec, level="info"):
    """Log a compact, readable summary table for `spec` using module logger."""
    table = inputs_markdown(inputs, spec=spec)
    logger.info("\nInputs Summary\n%s", table)


def save_inputs(inputs, path):
    """Save inputs to JSON for reproducibility (spec-agnostic)."""
    p = Path(path)
    if p.parent:
        p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(inputs, f, indent=2, sort_keys=True)


class InputsRegistry:
    """Bind a SPEC once and reuse methods (no globals required).

    Example:
        reg = InputsRegistry(SPEC_CALIBRATE)
        inputs = reg.collect(dataset_manifest, args, config_dict)
        reg.validate(inputs)
        reg.log(inputs)
        reg.save(inputs, Path(out) / "inputs.json")
    """

    def __init__(self, spec):
        self.spec = spec

    def collect(self, dataset_manifest, args, config_dict):
        return collect_inputs(dataset_manifest, args, config_dict, spec=self.spec)

    def validate(self, inputs):
        return validate_inputs(inputs, spec=self.spec)

    def markdown(self, inputs):
        return inputs_markdown(inputs, spec=self.spec)

    def log(self, inputs, level="info"):
        return log_inputs_summary(inputs, spec=self.spec, level=level)

    def save(self, inputs, path):
        return save_inputs(inputs, path)


# ---------- 2) Reuse helpers ----------

def compose_specs(*specs):
    """Merge multiple SPEC lists into one (later entries override by path)."""
    by_path = {}
    for sp in specs:
        for f in sp:
            path = f["path"]
            by_path[path] = {**by_path.get(path, {}), **f}
    return list(by_path.values())


def pick(inputs, keys):
    """Return a flat dict of selected dotted keys -> values from `inputs`."""
    out = {}
    for k in keys:
        out[k] = _nested_get(inputs, k)
    return out


# ---------- 3) Small utilities ----------

def _nested_set(d, dotted, value):
    parts = dotted.split(".")
    cur = d
    for p in parts[:-1]:
        if p not in cur or not isinstance(cur[p], dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value


def _nested_get(d, dotted, default=None):
    parts = dotted.split(".")
    cur = d
    for p in parts:
        if not isinstance(cur, dict) or p not in cur:
            return default
        cur = cur[p]
    return cur


def _is_missing(v):
    return v is None or (isinstance(v, str) and v.strip() == "")


def _is_number(v):
    try:
        float(v)
        return True
    except Exception:
        return False


def _is_int(v):
    try:
        return int(v) == v if isinstance(v, int) else float(v).is_integer()
    except Exception:
        return False


def _call_validator(fn, value):
    try:
        return bool(fn(value))
    except Exception:
        return False


def _safe_getattr(obj, name, default=None):
    try:
        return getattr(obj, name)
    except Exception:
        return default


def _fmt_value(v, maxlen=60):
    if isinstance(v, (list, tuple)):
        s = ", ".join([_short(repr(x), maxlen=20) for x in list(v)[:5]])
        if len(v) > 5:
            s += ", …"
        return "[" + s + "]"
    return _short(repr(v), maxlen=maxlen)


def _short(s, maxlen=60):
    s = str(s)
    return s if len(s) <= maxlen else s[: maxlen - 1] + "…"

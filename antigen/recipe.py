"""
recipe.py — lightweight "scripts as recipes" framework for antigen

Goals:
  • Make each script a *recipe*: Context (inputs), Operations (steps), Outputs.
  • Be clear, editable, and dependency-light.
  • Strong logging: human logs + JSONL run log, plus auto docs.

Key pieces:
  - Operation: single-responsibility step with declared needs/provides.
  - Recipe: ordered list of Operations + context spec (from antigen.inputs).
  - RunLogger: structured JSONL + standard logging integration.

"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import inspect
import json
import time
import logging
from typing import Dict, Any, List, Type

from antigen.fiber import load_fiber_positions
from antigen.wavelength import get_rectified_wavelength
from antigen.io import load_reduced_data, load_calspec_spectrum, read_extinction_table
from antigen.extinction import apply_extinction_correction
from antigen.calibrate import measure_response, build_psf_and_dar, extract_optimal_spectrum
from antigen.inputs import InputsRegistry, inputs_markdown
from antigen.plot import plot_spectrum_with_standard

logger = logging.getLogger("antigen.recipe")


# -----------------------------
# Operation base class
# -----------------------------
@dataclass
class Operation:
    """Base class for all operations in a recipe."""
    _name: str = field(init=False)
    _needs: List[str] = field(init=False, default_factory=list)
    _provides: List[str] = field(init=False, default_factory=list)
    _summary: str = field(init=False, default="")
    _init_params: List[str] = field(init=False, default_factory=list)

    def __post_init__(self):
        # Get initialization parameters (excluding self and optional params)
        sig = inspect.signature(self.__class__.__init__)
        self._init_params = [
            param.name for param in sig.parameters.values()
            if param.name != 'self' and param.default == param.empty
        ]

    def get_all_dependencies(self) -> List[str]:
        """Get both state requirements and initialization parameters."""
        return list(set(self._needs + self._init_params))

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        raise NotImplementedError



# -----------------------------
# Run logger (JSONL + std logging)
# -----------------------------
class RunLogger:
    def __init__(self, outdir: Path, recipe_name: str, *, pretty: bool = True, use_color: bool | None = None):
        self.outdir = outdir
        self.recipe_name = recipe_name
        self.jsonl_path = outdir / f"{recipe_name}.run.jsonl"
        self._fh = None
        self.pretty = pretty
        # auto-detect TTY for color unless explicitly set
        try:
            import sys
            self.use_color = use_color if use_color is not None else sys.stderr.isatty()
        except Exception:
            self.use_color = False
        self._ensure_file()

    def _ensure_file(self) -> None:
        self.outdir.mkdir(parents=True, exist_ok=True)
        self._fh = self.jsonl_path.open("a", encoding="utf-8")

    def emit(self, level: str, event: str, **kv: Any) -> None:
        # 1) Structured JSONL (unchanged)
        rec = {"ts": time.time(), "level": level, "event": event, **kv}
        self._fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
        self._fh.flush()

        # 2) Human log line
        lvl_no = getattr(logging, level.upper(), logging.INFO)
        if self.pretty:
            msg = self._pretty_line(event, kv)
        else:
            msg = f"{event} | " + " ".join(f"{k}={v!r}" for k, v in kv.items())
        logger.log(lvl_no, msg)

    def debug(self, event: str, **kv: Any) -> None: self.emit("debug", event, **kv)
    def info(self, event: str, **kv: Any) -> None:  self.emit("info", event, **kv)
    def warning(self, event: str, **kv: Any) -> None: self.emit("warning", event, **kv)
    def error(self, event: str, **kv: Any) -> None:   self.emit("error", event, **kv)

    def close(self) -> None:
        try:
            if self._fh: self._fh.close()
        finally:
            self._fh = None

    # ---- Pretty helpers ----
    def _pretty_line(self, event: str, kv: dict) -> str:
        icon = {
            "recipe_start": "🍳",
            "recipe_done":  "✅",
            "op_start":     "▶",
            "op_done":      "✓",
            "op_missing_inputs": "⚠",
            "error":        "✖",
        }.get(event, "•")

        # lift out common fields if present
        step   = kv.get("step")
        total  = kv.get("total")
        op     = kv.get("op")
        secs   = kv.get("seconds")

        # base headline
        head = []
        head.append(icon)
        if event in ("op_start", "op_done"):
            if step is not None and total is not None:
                head.append(f"{int(step)}/{int(total)}")
            elif step is not None:
                head.append(f"{int(step)}")
            if op:
                head.append(str(op))
        else:
            head.append(event)

        # tail bits
        tail = []
        if secs is not None:
            try:
                tail.append(f"{float(secs):.3f}s")
            except Exception:
                tail.append(f"{secs}")

        # remaining kvs, excluding ones we already surfaced
        exclude = {"step", "total", "op", "seconds"}
        rest = " ".join(f"{k}={kv[k]!r}" for k in kv if k not in exclude)
        if rest:
            tail.append(rest)

        # assemble
        line = " ".join(head)
        if tail:
            line += " " + " ".join(tail)

        # optional colorize simple levels
        if self.use_color:
            return self._colorize(event, line)
        return line

    def _colorize(self, event: str, s: str) -> str:
        # simple mapping; adjust if you want more nuance
        if event in ("recipe_start", "op_start"):
            code = 36   # cyan
        elif event in ("recipe_done", "op_done"):
            code = 32   # green
        elif event in ("op_missing_inputs", "warning"):
            code = 33   # yellow
        elif event in ("error",):
            code = 31   # red
        else:
            code = 0
        return f"\x1b[{code}m{s}\x1b[0m" if code else s


# -----------------------------
# Recipe container
# -----------------------------
@dataclass
class Recipe:
    name: str
    spec: List[Dict[str, Any]]  # Inputs spec used with InputsRegistry
    steps: List[Operation]
    outputs: List[str] = field(default_factory=list)  # keys expected in state
    description: str = ""

    def plan(self) -> str:
        """Return a human-readable plan (no execution)."""
        lines: List[str] = [f"Recipe: {self.name}"]
        if self.description:
            lines.append(self.description)
        lines.append("\nSteps:")
        for i, op in enumerate(self.steps, 1):
            need = ", ".join(op._needs) or "–"
            init = ", ".join(op._init_params) or "–"
            prov = ", ".join(op._provides) or "–"
            lines.append(
                f"  {i:02d}. {op._name}"
                f"\n      state needs: [{need}]"
                f"\n      init params: [{init}]"
                f"\n      provides: [{prov}]"
                f"\n      {op._summary}"
            )
        if self.outputs:
            lines.append("\nExpected outputs: " + ", ".join(self.outputs))
        return "\n".join(lines)

    def describe_markdown(self, manifest: Dict[str, Any], args, config_dict: Dict[str, Any]) -> str:
        """Produce Markdown docs: context table + step list."""
        reg = InputsRegistry(self.spec)
        inputs = reg.collect(manifest, args, config_dict)
        reg.validate(inputs)
        ctx_table = inputs_markdown(inputs, self.spec)
        lines: List[str] = [f"# Recipe: {self.name}"]
        if self.description:
            lines.append(self.description)
        lines.append("\n## Context (inputs)")
        lines.append(ctx_table)
        lines.append("\n## Steps")
        for i, op in enumerate(self.steps, 1):
            lines.append(f"{i}. **{op._name}** — {op._summary}")
        if self.outputs:
            lines.append("\n## Expected outputs")
            lines.append("- " + "\n- ".join(self.outputs))
        return "\n".join(lines)

    def run(self, manifest: Dict[str, Any], args, config_dict: Dict[str, Any], outdir: Path) -> Dict[str, Any]:
        """Execute the recipe end-to-end.

        Returns the final `state` dict containing outputs.
        """
        reg = InputsRegistry(self.spec)
        inputs = reg.collect(manifest, args, config_dict)
        reg.validate(inputs)
        reg.log(inputs)

        rlog = RunLogger(outdir, self.name)
        state: Dict[str, Any] = {"inputs": inputs, "outdir": Path(outdir)}

        rlog.info("recipe_start", name=self.name)
        for idx, op in enumerate(self.steps, 1):
            # Check needs
            missing = [k for k in op._needs if k not in state]
            if missing:
                rlog.error("op_missing_inputs", step=idx, op=op._name, missing=missing)
                raise ValueError(f"Operation {op._name} missing inputs: {missing}")
            rlog.info("op_start", step=idx, total=len(self.steps), op=op._name)
            t0 = time.time()
            op.run(state, rlog)
            dt = time.time() - t0
            rlog.info("op_done", step=idx, total=len(self.steps), op=op._name, seconds=round(dt, 3))
        rlog.info("recipe_done", name=self.name)
        rlog.close()
        return state


# -----------------------------
# Example operations (reusable)
# -----------------------------
@dataclass
class LoadFiberPositions(Operation):
    """Load fiber positions from instrument configuration."""
    unit_instrument: str
    unit_id: str
    ndithers: int
    dither_number: List[int]
    config_dict: Dict[str, Any]

    def __post_init__(self):
        super().__post_init__()
        self._name = "LoadFiberPositions"
        self._needs = []
        self._provides = ["fiber_x", "fiber_y"]
        self._summary = "Load fiber geometry for the instrument."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        fx, fy = load_fiber_positions(
            self.unit_instrument,
            self.ndithers,
            self.dither_number,
            self.config_dict
        )
        state.update({
            "fiber_x": fx,
            "fiber_y": fy
        })
        rlog.debug("fiber_positions_loaded")


@dataclass
class GetWavelengthGrid(Operation):
    """Generate wavelength grid from configuration."""
    config_dict: Dict[str, Any]

    def __post_init__(self):
        super().__post_init__()
        self._name = "GetWavelengthGrid"
        self._needs = []
        self._provides = ["def_wave"]
        self._summary = "Generate rectified wavelength grid."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        wave = get_rectified_wavelength(self.config_dict)
        state.update({
            "def_wave": wave
        })
        rlog.debug("wavelength_grid_created", n_wave=len(wave))


@dataclass
class LoadReducedData(Operation):
    """Load reduced spectra and error data."""
    in_folder: str
    reduced_files: List[str]

    def __post_init__(self):
        super().__post_init__()
        self._name = "LoadReducedData"
        self._needs = []
        self._provides = ["reduced_spectra", "reduced_error", "header"]
        self._summary = "Load reduced spectral frames and error data."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        rs, re, hdr = load_reduced_data(self.in_folder, self.reduced_files)
        state.update({
            "reduced_spectra": rs,
            "reduced_error": re,
            "header": hdr
        })
        rlog.debug("data_loaded", n_files=len(self.reduced_files))


@dataclass
class ModelPSFAndDAR(Operation):
    """Model PSF and DAR from spectral data."""
    extraction_radius: float
    fiber_radius: float

    def __post_init__(self):
        super().__post_init__()
        self._name = "ModelPSFAndDAR"
        self._needs = ["reduced_spectra", "reduced_error", "fiber_x", "fiber_y", "def_wave"]
        self._provides = ["modeling"]
        self._summary = "Fit PSF/DAR/FWHM as functions of wavelength."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:

        modeling = build_psf_and_dar(
            fiber_x=state["fiber_x"],
            fiber_y=state["fiber_y"],
            def_wave=state["def_wave"],
            reduced_spectra=state["reduced_spectra"],
            reduced_error=state["reduced_error"],
            extraction_radius=self.extraction_radius,
            fiber_radius=self.fiber_radius,
        )
        state["modeling"] = modeling
        rlog.debug("modeling_ready", have=list(modeling.keys()))


@dataclass
class ExtractSpectrum(Operation):
    """Extract optimal 1D spectrum using PSF-weighted fiber fluxes."""

    def __post_init__(self):
        super().__post_init__()
        self._name = "ExtractSpectrum"
        self._needs = [
            "reduced_spectra",
            "reduced_error",
            "modeling",  # Contains dar_model, sources, X, Y, measured_fwhm
            "fiber_x",
            "fiber_y",
            "def_wave"
        ]
        self._provides = ["spectrum", "spectrum_error"]
        self._summary = "Extract optimal 1D spectrum using PSF weights."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:

        modeling = state["modeling"]
        spectrum, spectrum_error = extract_optimal_spectrum(
            reduced_spectra=state["reduced_spectra"],
            reduced_error=state["reduced_error"],
            dar_model=modeling["dar_model"],
            sources=modeling["sources"],
            X=modeling["X"],
            Y=modeling["Y"],
            measured_fwhm=modeling["measured_fwhm"],
            fiber_x=state["fiber_x"],
            fiber_y=state["fiber_y"],
            psf_interp=modeling["psf_interp"],
            def_wave=state["def_wave"]
        )

        state.update({
            "spectrum": spectrum,
            "spectrum_error": spectrum_error
        })
        rlog.debug("spectrum_extracted", n_wave=len(spectrum))

@dataclass
class LoadExtinctionTable(Operation):
    """Load atmospheric extinction table."""
    
    def __post_init__(self):
        super().__post_init__()
        self._name = "LoadExtinctionTable"
        self._needs = []
        self._provides = ["extinction_table"]
        self._summary = "Load McDonald Observatory extinction curve."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        extinction_table = read_extinction_table()
        state["extinction_table"] = extinction_table
        rlog.debug("extinction_table_loaded")


@dataclass
class GetAirmass(Operation):
    """Extract airmass from header."""
    
    def __post_init__(self):
        super().__post_init__()
        self._name = "GetAirmass"
        self._needs = ["header"]
        self._provides = ["airmass"]
        self._summary = "Get airmass value from observation header."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        airmass = float(state["header"]["AIRMASS"])
        state["airmass"] = airmass
        rlog.debug("airmass_extracted", airmass=airmass)


@dataclass
class ApplyExtinctionCorrection(Operation):
    """Apply extinction correction to spectrum."""
    
    def __post_init__(self):
        super().__post_init__()
        self._name = "ApplyExtinctionCorrection"
        self._needs = [
            "def_wave",
            "spectrum",
            "spectrum_error",
            "extinction_table",
            "airmass"
        ]
        self._provides = ["spec_corr", "err_corr"]
        self._summary = "Apply atmospheric extinction correction to spectrum."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        spec_corr, err_corr = apply_extinction_correction(
            obs_wave=state["def_wave"],
            flux=state["spectrum"],
            flux_err=state["spectrum_error"],
            ext_wave=state["extinction_table"]["wavelength"],
            A_lambda_airmass=state["extinction_table"]["mags_per_airmass"],
            airmass=state["airmass"]
        )
        state.update({
            "spec_corr": spec_corr,
            "err_corr": err_corr
        })
        rlog.info("extinction_correction_applied", airmass=state["airmass"])


@dataclass
class LoadCalibrationSpectrum(Operation):
    """Load standard star calibration spectrum."""
    standard_name: str

    def __post_init__(self):
        super().__post_init__()
        self._name = "LoadCalibrationSpectrum"
        self._needs = []
        self._provides = ["cal_spectrum"]
        self._summary = "Load CALSPEC standard star spectrum."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        cal_spectrum = load_calspec_spectrum(self.standard_name)
        state["cal_spectrum"] = cal_spectrum
        rlog.debug("calibration_spectrum_loaded", standard=self.standard_name)


@dataclass
class SaveResponsePlot(Operation):
    """Generate and save response function plot."""
    output_folder: Path

    def __post_init__(self):
        super().__post_init__()
        self._name = "SaveResponsePlot"
        self._needs = [
            "def_wave",
            "spec_corr",
            "err_corr",
            "cal_spectrum",
            "response"
        ]
        self._provides = ["flux_cal"]
        self._summary = "Save diagnostic plot of response function."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        plot_spectrum_with_standard(
            state["flux_calibrated"],
            state["error_calibrated"],
            state["def_wave"],
            state["cal_spectrum"]["WAVELENGTH"],
            state["cal_spectrum"]["FLUX"],
            1.0 / state["response"],
            outfolder=self.output_folder
        )
        state["flux_cal"] = True
        rlog.info("response_plot_saved", output_folder=str(self.output_folder))

@dataclass
class MeasureResponse(Operation):
    """Measure instrument response function."""
    
    def __post_init__(self):
        super().__post_init__()
        self._name = "MeasureResponse"
        self._needs = ["def_wave", "spec_corr", "cal_spectrum"]
        self._provides = ["response"]
        self._summary = "Compute instrument response from standard star observation."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        response = measure_response(
            obs_wave=state["def_wave"],
            obs_flux=state["spec_corr"],
            std_wave=state["cal_spectrum"]["WAVELENGTH"],
            std_flux=state["cal_spectrum"]["FLUX"]
        )
        state["response"] = response
        rlog.debug("response_measured")


@dataclass
class ApplyResponse(Operation):
    """Apply response correction to spectrum."""
    
    def __post_init__(self):
        super().__post_init__()
        self._name = "ApplyResponse"
        self._needs = ["spec_corr", "err_corr", "response"]
        self._provides = ["flux_calibrated", "error_calibrated"]
        self._summary = "Apply response correction to spectrum and error."

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        spectrum = state["spec_corr"] * state["response"]
        error = state["err_corr"] * state["response"]
        
        state["flux_calibrated"] = spectrum
        state["error_calibrated"] = error
        rlog.debug("response_applied")
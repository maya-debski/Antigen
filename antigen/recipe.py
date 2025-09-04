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
import json
import time
import logging
from typing import Dict, Any, List, Type

from antigen.inputs import InputsRegistry, inputs_markdown  # reuse your spec tools
from antigen import config

logger = logging.getLogger("antigen.recipe")


# -----------------------------
# Operation base class
# -----------------------------
@dataclass
class Operation:
    name: str
    needs: List[str] = field(default_factory=list)     # keys required in state
    provides: List[str] = field(default_factory=list)  # keys added to state
    summary: str = ""

    def run(self, state: Dict[str, Any], rlog: "RunLogger") -> None:
        """Implement in subclass: mutate `state` with provided keys.
        Should log progress via `rlog.info()` etc.
        """
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
            need = ", ".join(op.needs) or "–"
            prov = ", ".join(op.provides) or "–"
            lines.append(f"  {i:02d}. {op.name} | needs: [{need}] -> provides: [{prov}]  {op.summary}")
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
            lines.append(f"{i}. **{op.name}** — {op.summary}")
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
            missing = [k for k in op.needs if k not in state]
            if missing:
                rlog.error("op_missing_inputs", step=idx, op=op.name, missing=missing)
                raise ValueError(f"Operation {op.name} missing inputs: {missing}")
            rlog.info("op_start", step=idx, total=len(self.steps), op=op.name)
            t0 = time.time()
            op.run(state, rlog)
            dt = time.time() - t0
            rlog.info("op_done", step=idx, total=len(self.steps), op=op.name, seconds=round(dt, 3))
        rlog.info("recipe_done", name=self.name)
        rlog.close()
        return state


# -----------------------------
# Example operations (reusable)
# -----------------------------
class PrepareData(Operation):
    def __init__(self):
        super().__init__(
            name="PrepareData",
            needs=["ctx"],
            provides=["fiber_x", "fiber_y", "def_wave", "reduced_spectra", "reduced_error", "header"],
            summary="Load fiber geometry, rectified wavelength, reduced frames.",
        )

    def run(self, state: Dict[str, Any], rlog: RunLogger) -> None:
        from antigen import fiber, wavelength, io
        ctx = state["ctx"]
        fx, fy = fiber.load_fiber_positions(ctx.unit_instrument, ctx.ndithers, ctx.dither_number, ctx.config_dict)
        wave = wavelength.get_rectified_wavelength(ctx.config_dict)
        rs, re, hdr = io.load_reduced_data(ctx.in_folder, ctx.reduced_files)
        state.update({"fiber_x": fx, "fiber_y": fy, "def_wave": wave, "reduced_spectra": rs,
                      "reduced_error": re, "header": hdr})
        rlog.debug("data_prepared", n_wave=len(wave), n_files=len(ctx.reduced_files))


class ModelPSFAndDAR(Operation):
    def __init__(self):
        super().__init__(
            name="ModelPSFAndDAR",
            needs=["ctx", "reduced_spectra", "reduced_error", "fiber_x", "fiber_y", "def_wave"],
            provides=["modeling"],
            summary="Fit PSF/DAR/FWHM as functions of wavelength.",
        )

    def run(self, state: Dict[str, Any], rlog: RunLogger) -> None:
        from antigen.calibrate import build_psf_and_dar
        ctx = state["ctx"]

        modeling = build_psf_and_dar(
            fiber_x=state["fiber_x"],
            fiber_y=state["fiber_y"],
            def_wave=state["def_wave"],
            reduced_spectra=state["reduced_spectra"],
            reduced_error=state["reduced_error"],
            extraction_radius=ctx.extraction_radius,
            fiber_radius=ctx.fiber_radius,
        )
        state["modeling"] = modeling
        rlog.debug("modeling_ready", have=list(modeling.keys()))


class ExtractSpectrum(Operation):
    def __init__(self):
        super().__init__(
            name="ExtractSpectrum",
            needs=["ctx", "reduced_spectra", "reduced_error", "modeling", "fiber_x", "fiber_y", "def_wave"],
            provides=["spectrum", "spectrum_error"],
            summary="Optimal extraction along the DAR track.",
        )

    def run(self, state: Dict[str, Any], rlog: RunLogger) -> None:
        from antigen.calibrate import extract_optimal_spectrum
        spec, err = extract_optimal_spectrum(
            reduced_spectra=state["reduced_spectra"],
            reduced_error=state["reduced_error"],
            dar_model=state["modeling"]["dar_model"],
            sources=state["modeling"].get("sources"),
            X=state["modeling"].get("X"),
            Y=state["modeling"].get("Y"),
            measured_fwhm=state["modeling"].get("measured_fwhm"),
            fiber_x=state["fiber_x"],
            fiber_y=state["fiber_y"],
            psf_interp=state["modeling"].get("psf_interp"),
            def_wave=state["def_wave"]
        )
        state["spectrum"] = spec
        state["spectrum_error"] = err
        rlog.debug("spectrum_extracted", n=len(spec))


class ApplyExtinction(Operation):
    def __init__(self):
        super().__init__(
            name="ApplyExtinction",
            needs=["ctx", "spectrum", "spectrum_error", "def_wave", "header"],
            provides=["spec_corr", "err_corr"],
            summary="Correct for site extinction using header AIRMASS.",
        )

    def run(self, state: Dict[str, Any], rlog: RunLogger) -> None:
        from antigen import extinction, io
        xtab = io.read_extinction_table()
        airmass = float(state["header"].get("AIRMASS", 1.0))
        factor = extinction.compute_extinction_factor(state["def_wave"], xtab["wavelength"],
                                                     xtab["mags_per_airmass"], airmass)
        state["spec_corr"] = state["spectrum"] * factor
        state["err_corr"] = state["spectrum_error"] * factor
        rlog.debug("extinction_applied", airmass=airmass)


class MeasureResponse(Operation):
    def __init__(self):
        super().__init__(
            name="MeasureResponse",
            needs=["ctx", "def_wave", "spec_corr", "err_corr"],
            provides=["flux_cal", "flux_cal_err", "response"],
            summary="Measure response vs. standard star and apply calibration.",
        )


    def run(self, state: Dict[str, Any], rlog: RunLogger) -> None:
        from antigen.calibrate import measure_and_apply_response
        from antigen import io
        ctx = state["ctx"]
        cal_star = io.load_calspec_spectrum(ctx.standard_name)
        flux_cal, flux_cal_err, response = measure_and_apply_response(
            state["def_wave"], state["spec_corr"], state["err_corr"], cal_star, ctx.output_folder, window=251
        )
        state["flux_cal"] = flux_cal
        state["flux_cal_err"] = flux_cal_err
        state["response"] = response
        # io.save_response_curve(ctx.output_folder, state["def_wave"], response)
        rlog.debug("response_measured", out=str(ctx.output_folder))

@dataclass
class BuildContext:
    """Base context class that provides common context building functionality."""
    inputs: Dict[str, Any]
    
    def __post_init__(self):
        self._process_inputs()
    
    def _process_inputs(self):
        """Process inputs according to their categories.
        Override this in subclasses to add specific processing."""
        pass

    def ensure_paths(self):
        """Ensure all path-like attributes are properly converted to Path objects"""
        for attr_name, attr_value in self.__dict__.items():
            if isinstance(attr_value, str) and ('path' in attr_name.lower() or 
                                              'folder' in attr_name.lower() or 
                                              'dir' in attr_name.lower()):
                setattr(self, attr_name, Path(attr_value).expanduser())

class CalibrateContext(BuildContext):
    """Specific context for calibration operations."""
    
    def _process_inputs(self):
        # CLI inputs
        self.standard_name = self.inputs["cli"]["standard_name"]
        self.output_folder = Path(self.inputs["cli"]["output_folder"]).expanduser()
        self.extraction_radius = float(self.inputs["cli"].get("extraction_radius", 2.0))

        # Dataset inputs
        self.unit_instrument = self.inputs["dataset"]["unit_instrument"]
        self.unit_id = self.inputs["dataset"]["unit_id"]
        self.ndithers = int(self.inputs["dataset"]["ndithers"])
        self.dither_number = [int(n) for n in self.inputs["dataset"]["dither_number"]]
        self.in_folder = Path(self.inputs["dataset"]["in_folder"]).expanduser()
        self.reduced_files = [str(self.in_folder / f) for f in self.inputs["dataset"]["reduced_files"]]

        # Config inputs
        self.fiber_radius = float(self.inputs["config"]["fiber_radius"])
        
        # Construct config dict
        self.config_dict = config.build_config_for_element(
            self.unit_instrument.lower(), self.unit_id.upper()
        )
        
        # Ensure output folder exists
        self.output_folder.mkdir(parents=True, exist_ok=True)

class BuildGenericContext(Operation):
    """Generic context building operation that can work with any context type."""
    
    def __init__(self, context_class: Type[BuildContext]):
        super().__init__(
            name=f"Build{context_class.__name__}",
            needs=["inputs"],
            provides=["ctx"],
            summary=f"Create {context_class.__name__} from validated inputs.",
        )
        self.context_class = context_class

    def run(self, state: Dict[str, Any], rlog: RunLogger) -> None:
        ctx = self.context_class(state["inputs"])
        state["ctx"] = ctx
        rlog.debug("ctx_ready", context_type=self.context_class.__name__)
"""Self-contained HTML rendering of one sealed shadow record.

The report is a rendering of a record, never a recomputation of one: every
number shown is read from the sealed JSON (and, when supplied, its validation
JSON). Output is a single file with no external requests, so a report stays
readable from an archive with no network and no installed viewer.
"""

from __future__ import annotations

from html import escape
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import __version__


CSS = """
:root{
  --bg:#f3f5f7; --surface:#ffffff; --surface-2:#eaeef2; --border:#d7dee5;
  --ink:#16202b; --ink-muted:#5c6b7a; --ink-faint:#8b98a6;
  --accent:#b9791c; --accent-ink:#7a4f0f;
  --good:#0f7a52; --good-bg:#e3f5ee;
  --critical:#c5433f; --critical-bg:#fbe9e8;
  --trace:#3a5da8; --trace-bg:#e8edf7;
  --radius:10px; color-scheme:light;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#0a0e13; --surface:#10161f; --surface-2:#182231; --border:#26313f;
    --ink:#e8edf2; --ink-muted:#8ea0b5; --ink-faint:#5f7185;
    --accent:#d9a441; --accent-ink:#f0c374;
    --good:#3ecf8e; --good-bg:#123528;
    --critical:#ef6b67; --critical-bg:#3a1a1a;
    --trace:#7ea1e8; --trace-bg:#182338;
    color-scheme:dark;
  }
}
:root[data-theme="dark"]{
  --bg:#0a0e13; --surface:#10161f; --surface-2:#182231; --border:#26313f;
  --ink:#e8edf2; --ink-muted:#8ea0b5; --ink-faint:#5f7185;
  --accent:#d9a441; --accent-ink:#f0c374;
  --good:#3ecf8e; --good-bg:#123528;
  --critical:#ef6b67; --critical-bg:#3a1a1a;
  --trace:#7ea1e8; --trace-bg:#182338;
  color-scheme:dark;
}
*{box-sizing:border-box;}
body{margin:0; background:var(--bg); color:var(--ink); line-height:1.45;
  font-family:ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
.mono{font-family:ui-monospace,"SF Mono","Cascadia Code",Consolas,monospace;
  font-variant-numeric:tabular-nums;}
.wrap{max-width:960px; margin:0 auto; padding:28px 24px 56px;}
.topbar{display:flex; justify-content:space-between; align-items:baseline; gap:20px;
  flex-wrap:wrap; padding-bottom:18px; border-bottom:1px solid var(--border);
  margin-bottom:22px;}
.mark{font-family:Georgia,"Iowan Old Style",ui-serif,serif; font-weight:700;
  font-size:1.5rem;}
.mark .dim{color:var(--ink-faint); font-weight:400;}
.sub{font-size:.72rem; text-transform:uppercase; letter-spacing:.12em;
  color:var(--ink-muted);}
.badge{display:inline-flex; align-items:center; gap:7px; padding:5px 11px;
  border-radius:999px; border:1px solid var(--border); background:var(--surface-2);
  font-size:.74rem; text-transform:uppercase; letter-spacing:.08em; font-weight:600;}
.badge.good{color:var(--good); background:var(--good-bg); border-color:var(--good);}
.badge.critical{color:var(--critical); background:var(--critical-bg);
  border-color:var(--critical);}
.badge.accent{color:var(--accent-ink); border-color:var(--accent);}
.card{background:var(--surface); border:1px solid var(--border);
  border-radius:var(--radius); padding:18px 20px; margin-bottom:18px;}
.card > h2{margin:0 0 4px; font-size:.78rem; text-transform:uppercase;
  letter-spacing:.1em; color:var(--ink-muted);}
.card > .note{margin:0 0 14px; font-size:.82rem; color:var(--ink-faint);}
.grid{display:grid; grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); gap:14px;}
.stat .k{font-size:.7rem; text-transform:uppercase; letter-spacing:.08em;
  color:var(--ink-muted);}
.stat .v{font-size:1.22rem; font-weight:650;}
.stat .v.small{font-size:.9rem; font-weight:500; word-break:break-all;}
table{width:100%; border-collapse:collapse; font-size:.88rem;}
td{padding:7px 0; border-bottom:1px solid var(--border); vertical-align:top;}
td:first-child{color:var(--ink-muted); padding-right:16px;}
td:last-child{text-align:right;}
tr:last-child td{border-bottom:none;}
.plot{overflow-x:auto;}
svg{display:block; max-width:100%; height:auto;}
.pass{color:var(--good); font-weight:600;}
.fail{color:var(--critical); font-weight:600;}
.na{color:var(--ink-faint);}
.scope{font-size:.82rem; color:var(--ink-muted); border-left:3px solid var(--accent);
  padding:8px 0 8px 14px; margin:0;}
footer{margin-top:22px; padding-top:14px; border-top:1px solid var(--border);
  font-size:.76rem; color:var(--ink-faint);}
"""


def _number(value: Any, digits: int = 3, suffix: str = "") -> str:
    if value is None:
        return "&mdash;"
    try:
        return f"{float(value):.{digits}f}{suffix}"
    except (TypeError, ValueError):
        return escape(str(value))


def _text(value: Any, fallback: str = "&mdash;") -> str:
    if value is None or value == "":
        return fallback
    return escape(str(value))


def _short_hash(value: Any, keep: int = 16) -> str:
    if not isinstance(value, str) or not value:
        return "&mdash;"
    return escape(value[:keep]) + ("&hellip;" if len(value) > keep else "")


def _stat(key: str, value: str) -> str:
    return (
        '<div class="stat">'
        f'<div class="k">{key}</div>'
        f'<div class="v mono">{value}</div></div>'
    )


def _rows(pairs: Sequence[tuple[str, str]]) -> str:
    body = "".join(
        f"<tr><td>{key}</td><td class='mono'>{value}</td></tr>" for key, value in pairs
    )
    return f"<table>{body}</table>"


def _identity_line(record: Mapping[str, Any]) -> str:
    identity = record.get("run_identity") or {}
    protein = identity.get("protein")
    pocket = identity.get("pocket")
    if protein and pocket:
        return f"{escape(str(protein))} : {escape(str(pocket))}"
    return escape(str(record.get("run_id", "unknown run")))


def render_trace_svg(
    series: Mapping[str, Any], threshold: float | None, feature_mean: float | None
) -> str:
    """Draw the prefix RMSD trace as inline SVG.

    Inline because a report must stay self-contained: no plotting library at
    render time and no image file to lose beside the HTML.
    """
    times = list(series.get("time_ns") or ())
    rmsd = list(series.get("corrected_pose_rmsd_angstrom") or ())
    window = list(series.get("in_feature_window") or ())
    if len(times) < 2 or len(rmsd) != len(times):
        return '<p class="note">This record predates per-frame series capture.</p>'

    width, height = 640, 240
    left, right, top, bottom = 56, 16, 16, 34
    x_max = max(times)
    y_max = max([*rmsd, threshold or 0.0]) * 1.15 or 1.0

    def x_of(value: float) -> float:
        return left + (value / x_max) * (width - left - right)

    def y_of(value: float) -> float:
        return height - bottom - (value / y_max) * (height - top - bottom)

    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        'role="img" aria-label="Corrected ligand pose RMSD against time">'
    ]
    if window and len(window) == len(times) and any(window):
        inside = [time for time, flag in zip(times, window) if flag]
        x0, x1 = x_of(min(inside)), x_of(max(inside))
        parts.append(
            f'<rect x="{x0:.1f}" y="{top}" width="{max(x1 - x0, 1):.1f}" '
            f'height="{height - top - bottom:.1f}" fill="var(--trace-bg)"/>'
        )
    for step in range(5):
        value = y_max * step / 4
        y = y_of(value)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width - right}" y2="{y:.1f}" '
            'stroke="var(--border)" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" '
            f'font-size="11" fill="var(--ink-faint)">{value:.1f}</text>'
        )
    for step in range(5):
        value = x_max * step / 4
        x = x_of(value)
        parts.append(
            f'<text x="{x:.1f}" y="{height - bottom + 18:.1f}" text-anchor="middle" '
            f'font-size="11" fill="var(--ink-faint)">{value:g}</text>'
        )
    parts.append(
        f'<text x="{width - right}" y="{height - 4}" text-anchor="end" '
        'font-size="11" fill="var(--ink-faint)">ns</text>'
    )
    if threshold is not None and threshold <= y_max:
        y = y_of(threshold)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width - right}" y2="{y:.1f}" '
            'stroke="var(--critical)" stroke-width="1.5" stroke-dasharray="6 4"/>'
        )
        parts.append(
            f'<text x="{left + 6}" y="{y - 6:.1f}" font-size="11" '
            f'fill="var(--critical)">threshold {threshold:.4f} &#8491;</text>'
        )
    if feature_mean is not None and feature_mean <= y_max:
        y = y_of(feature_mean)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width - right}" y2="{y:.1f}" '
            'stroke="var(--accent)" stroke-width="1.5"/>'
        )
    points = " ".join(
        f"{x_of(time):.1f},{y_of(value):.1f}" for time, value in zip(times, rmsd)
    )
    parts.append(
        f'<polyline fill="none" stroke="var(--trace)" stroke-width="2.2" '
        f'stroke-linejoin="round" stroke-linecap="round" points="{points}"/>'
    )
    parts.append("</svg>")
    return f'<div class="plot">{"".join(parts)}</div>'


def _decision_card(record: Mapping[str, Any]) -> str:
    policy_type = record.get("policy_type", "standardized_logistic")
    measurement = record.get("measurement") or {}
    feature_mean = measurement.get("corrected_pose_rmsd_mean_angstrom")
    stats = [
        _stat(
            "Corrected pose RMSD", _number(feature_mean, 4, " &#8491;")
        ),
        _stat("Forecast", _text(record.get("forecast"))),
        _stat("Recommendation", _text(record.get("recommendation"))),
    ]
    if policy_type == "threshold_rule":
        margin = record.get("signed_stop_margin_angstrom")
        sign = "+" if isinstance(margin, (int, float)) and margin > 0 else ""
        stats.insert(
            1,
            _stat(
                "Stop threshold",
                _number(record.get("stop_threshold_angstrom"), 4, " &#8491;"),
            ),
        )
        stats.insert(2, _stat("Signed margin", sign + _number(margin, 4, " &#8491;")))
    else:
        stats.insert(
            1,
            _stat(
                "Uncalibrated score",
                _number(record.get("uncalibrated_retention_score"), 6),
            ),
        )
        stats.insert(
            2,
            _stat(
                "Decision threshold",
                _number(record.get("decision_threshold"), 3),
            ),
        )
    note = (
        "One-sided early-stop screen. Not meeting the stop condition is not a "
        "prediction that the pose is retained."
        if policy_type == "threshold_rule"
        else "The score is an uncalibrated ranking score, not a probability of "
        "binding or of drug activity."
    )
    window = (
        f"({_number(measurement.get('window_start_ns'), 1)}, "
        f"{_number(measurement.get('window_end_ns'), 1)}] ns"
    )
    return (
        '<section class="card"><h2>Frozen decision</h2>'
        f'<p class="note">{escape(str(record.get("policy_id", "")))} '
        f"&middot; {policy_type} &middot; feature window {window}</p>"
        f'<div class="grid">{"".join(stats)}</div>'
        f'<p class="scope" style="margin-top:16px">{note}</p></section>'
    )


def _quality_card(record: Mapping[str, Any]) -> str:
    quality = record.get("quality_control")
    if isinstance(quality, str):  # schema 1.0 records carried a bare verdict
        return (
            '<section class="card"><h2>Applicability &amp; QC</h2>'
            f'<p class="note">Schema 1.0 record: a single verdict, no named checks.</p>'
            f"{_rows([('Quality control', escape(quality))])}</section>"
        )
    if not isinstance(quality, Mapping):
        return ""
    checks = quality.get("checks") or {}
    labels = {
        "periodic_box_vectors_valid": "Periodic box vectors valid",
        "protein_ligand_selections_disjoint": "Protein / ligand selections",
        "checkpoint_coverage_complete": "Checkpoint coverage",
        "feature_window_populated": "Feature window populated",
        "ligand_within_declared_domain": "Domain profile",
        "future_frames_accessed": "Future frames accessed",
    }
    pairs = []
    for key, label in labels.items():
        if key not in checks:
            continue
        value = checks[key]
        if value is None:
            pairs.append((label, '<span class="na">not configured</span>'))
        elif key == "future_frames_accessed":
            css, word = ("pass", "NO") if value is False else ("fail", "YES")
            pairs.append((label, f'<span class="{css}">{word}</span>'))
        else:
            css, word = ("pass", "PASS") if value else ("fail", "FAIL")
            pairs.append((label, f'<span class="{css}">{word}</span>'))
    domain = quality.get("domain_profile") or {}
    fraction = domain.get("ligand_diameter_box_fraction")
    limit = domain.get("max_ligand_diameter_box_fraction")
    if fraction is not None:
        detail = _number(fraction, 4)
        if limit is not None:
            detail += f" (limit {_number(limit, 2)})"
        pairs.append(("Ligand diameter / box", detail))
    return (
        '<section class="card"><h2>Applicability &amp; QC</h2>'
        '<p class="note">Verified before the policy was applied &mdash; not after.</p>'
        f"{_rows(pairs)}</section>"
    )


def _record_card(record: Mapping[str, Any]) -> str:
    software = record.get("software") or {}
    pairs = [
        ("Forecast sealed", _text(record.get("created_utc"))),
        ("Record id", _text(record.get("record_id"))),
        (
            "Configuration hash",
            _short_hash(record.get("configuration_scientific_sha256")),
        ),
        (
            "Trajectory-prefix hash",
            _short_hash(record.get("trajectory_prefix_coordinates_sha256")),
        ),
        ("Topology hash", _short_hash(record.get("topology_sha256"))),
        (
            "Software",
            _text("posegate-md " + str(software.get("version", ""))),
        ),
    ]
    return (
        '<section class="card"><h2>Immutable record</h2>'
        '<p class="note">Sealed at capture; never edited in place.</p>'
        f"{_rows(pairs)}</section>"
    )


def _outcome_card(validation: Mapping[str, Any] | None) -> str:
    if validation is None:
        return (
            '<section class="card"><h2>Finalized outcome</h2>'
            '<p class="note">No validation record supplied. The late window has not '
            "been compared with this forecast in this report.</p></section>"
        )
    outcome = validation.get("outcome") or {}
    integrity = validation.get("record_integrity") or {}
    retained = outcome.get("pose_retained")
    category = str(validation.get("audit_category", ""))
    correct = validation.get("decision_was_correct")
    if correct is True:
        badge = f'<span class="badge good">{escape(category)}</span>'
    elif correct is False:
        badge = f'<span class="badge critical">{escape(category)}</span>'
    else:
        badge = f'<span class="badge">{escape(category)}</span>'
    intact = integrity.get("sealed_forecast_still_follows_from_sealed_measurement")
    if intact is None:
        modified = '<span class="na">&mdash;</span>'
    else:
        modified = (
            '<span class="pass">NO</span>'
            if intact
            else '<span class="fail">YES</span>'
        )
    pairs = [
        (
            "Late median RMSD",
            _number(outcome.get("late_pose_rmsd_median_angstrom"), 3, " &#8491;"),
        ),
        (
            "Retention boundary",
            _number(outcome.get("retained_rmsd_threshold_angstrom"), 3, " &#8491;"),
        ),
        (
            "Late window",
            f"({_number(outcome.get('late_window_start_ns'), 0)}, "
            f"{_number(outcome.get('late_window_end_ns'), 0)}] ns",
        ),
        ("Outcome", "Retained" if retained else "Non-retained"),
        ("Audit category", badge),
        ("Original forecast modified", modified),
    ]
    window = validation.get("feature_window_cross_check") or {}
    if window and window.get("matches_record") is False:
        pairs.append(
            (
                "Feature window on revalidation",
                '<span class="fail">'
                f"{window.get('frames_in_record')} &rarr; "
                f"{window.get('frames_on_revalidation')} frames</span>",
            )
        )
    elif window:
        pairs.append(
            (
                "Feature window on revalidation",
                f'<span class="pass">{window.get("frames_on_revalidation")} '
                "frames, unchanged</span>",
            )
        )
    read_only = integrity.get("record_file_read_only")
    if read_only is not None:
        pairs.append(
            (
                "Record file",
                '<span class="pass">read-only</span>'
                if read_only
                else '<span class="fail">writable</span>',
            )
        )
    return (
        '<section class="card"><h2>Finalized outcome &mdash; after the late window</h2>'
        f"{_rows(pairs)}</section>"
    )


def render_record_html(
    record: Mapping[str, Any], validation: Mapping[str, Any] | None = None
) -> str:
    identity = record.get("run_identity") or {}
    measurement = record.get("measurement") or {}
    series = measurement.get("series") or {}
    threshold = record.get("stop_threshold_angstrom")
    feature_mean = measurement.get("corrected_pose_rmsd_mean_angstrom")
    subtitle = " &middot; ".join(
        escape(str(identity[key]))
        for key in ("ligand", "replica")
        if identity.get(key)
    )
    status = (
        '<span class="badge accent">Outcome-blinded shadow</span>'
        if validation is None
        else '<span class="badge">Outcome finalized</span>'
    )
    checkpoint = _number(measurement.get("checkpoint_ns"), 1)
    plot = render_trace_svg(
        series,
        threshold if isinstance(threshold, (int, float)) else None,
        feature_mean if isinstance(feature_mean, (int, float)) else None,
    )
    return f"""<title>PoseGate-MD shadow record &mdash; {_identity_line(record)}</title>
<style>{CSS}</style>
<div class="wrap">
  <div class="topbar">
    <div>
      <div class="mark">PoseGate<span class="dim">-MD</span></div>
      <div class="sub">Shadow record &middot; {_identity_line(record)}</div>
    </div>
    <div>{status}</div>
  </div>

  <section class="card">
    <h2>Live checkpoint &mdash; {checkpoint} ns capture</h2>
    <p class="note">Corrected ligand-pose RMSD relative to the reference frame.
      Shaded band is the feature window the policy actually reads.
      {subtitle}</p>
    {plot}
  </section>

  {_decision_card(record)}
  {_quality_card(record)}
  {_record_card(record)}
  {_outcome_card(validation)}

  <p class="scope">{escape(str(record.get("scope_warning", "")))}</p>
  <footer>Rendered by posegate-md {escape(__version__)} from the sealed record.
    Every value on this page is read from the record;
    nothing is recomputed here.</footer>
</div>
"""


def write_record_report(
    output: str | Path,
    record: Mapping[str, Any],
    validation: Mapping[str, Any] | None = None,
) -> Path:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_record_html(record, validation), encoding="utf-8")
    return path


def load_validation(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict) or "outcome" not in value:
        raise ValueError(f"not a posegate validation record: {path}")
    return value

#!/usr/bin/env python3
r"""Plot four GNole correction ratios versus lambda = alpha_s * ln(Q/ET).

Example (Matplotlib is the only plotting dependency):
  python3 plot_corrections.py --hard hard.dat --single single.dat \
      --double double.dat --obs obs.dat --ll ll.dat --alphas 0.1 \
      --out corrections_lambda.png

Inputs are GNole output files, not raw differential histograms. Hard, single,
double and LL use the stored lambda cumulants; obs uses the already-normalized
lambda cumulant correction at its shifted bin centers. No integration,
extra coupling factor, flow-bin deposit or per-file renormalization is added.

The LL file is used exactly as supplied. With physical LL and SL component
runs, the first three ratios also contain a measurement difference. Use an
LL -sl-obs file if you want a common SL baseline. Obs is a finite-bin average,
not an exact pointwise correction. Per-file errors are not combined because
shared-seed ratios require cross-run covariance.

The script also saves the plotted values as CSV. --no-plot needs only the
Python standard library and writes just the CSV.
"""
import argparse
import csv
import math
from pathlib import Path
import re
import sys


def read_curve(path, correction=False):
    """Read a lambda cumulant or shifted observable-correction table."""
    path = Path(path)
    target = ("# observable correction (lambda)" if correction
              else "# dSdlambda (cumulative)")
    rows, metadata, active, found = [], {}, False, False
    text = path.read_text()
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        match = re.match(r"# ([^=]+?)\s*=\s*(\S+)", line)
        if match:
            metadata[match.group(1).strip()] = match.group(2)
        if line.startswith(("# dS", "# observable correction (")):
            active = line == target
            if active:
                if found:
                    raise ValueError(f"{path}: duplicate {target}")
                found = True
            continue
        if not active or not line or line.startswith("#"):
            continue
        try:
            values = list(map(float, line.split()))
        except ValueError:
            raise ValueError(f"{path}:{number}: invalid numeric row") from None
        expected = (4, 5) if correction else (2, 3)
        if len(values) not in expected or not all(map(math.isfinite, values)):
            raise ValueError(f"{path}:{number}: invalid or nonfinite table row")
        if len(values) == expected[1] and values[-1] < 0:
            raise ValueError(f"{path}:{number}: negative standard error")
        if correction and not values[0] < values[1] < values[2]:
            raise ValueError(f"{path}:{number}: invalid correction bin")
        rows.append((values[1], values[3]) if correction else tuple(values[:2]))
    if not rows:
        raise ValueError(f"{path}: missing or empty {target}")
    rows.sort()
    if any(b[0] <= a[0] for a, b in zip(rows, rows[1:])):
        raise ValueError(f"{path}: duplicate lambda coordinates")
    # The numeric header uses six significant digits. Prefer the exact -as
    # argument in the recorded command, when available.
    couplings = re.findall(r"(?:^|\s)-as\s+(\S+)", text.splitlines()[0])
    if couplings:
        exact = float(couplings[-1])
        if "alphas(rts)" in metadata and not math.isclose(
                exact, float(metadata["alphas(rts)"]), rel_tol=5e-6, abs_tol=1e-12):
            raise ValueError(f"{path}: recorded -as conflicts with alphas(rts)")
        metadata["alphas(rts)"] = str(exact)
    metadata["_sl_observable"] = ("-sl-obs" in text.split()
                                  or metadata.get("sl_observable") == "1")
    return {"x": [x for x, y in rows], "y": [y for x, y in rows],
            "metadata": metadata}


def build_ratios(curves, alphas):
    """Use LL edge coordinates; do not interpolate or integrate obs."""
    if not math.isfinite(alphas) or alphas <= 0:
        raise ValueError("alphas must be finite and positive")
    x, baseline = curves["ll"]["x"], curves["ll"]["y"]
    for name, curve in curves.items():
        if len(curve["x"]) != len(x) or any(
                not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-10)
                for a, b in zip(x, curve["x"])):
            raise ValueError(f"{name}: lambda grid does not match LL")
    ratios = {"lambda": x}
    for name in ("hard", "single", "double", "obs"):
        ratios[name] = [
            ((value - ll) if name != "obs" else value) / (ll * alphas)
            if ll > 0 else math.nan
            for value, ll in zip(curves[name]["y"], baseline)]
    if not any(ll > 0 for ll in baseline):
        raise ValueError("LL has no positive denominators")
    return ratios


def check_parameters(curves, alphas):
    """Different event counts/seeds are OK; physics settings must agree."""
    for key in ("alphas(rts)", "xmur", "xQ", "dy", "lnktmax", "etamax", "p"):
        known = [(name, float(curve["metadata"][key]))
                 for name, curve in curves.items() if key in curve["metadata"]]
        if known and any(not math.isclose(value, known[0][1], rel_tol=1e-8,
                                         abs_tol=1e-12)
                         for name, value in known):
            raise ValueError(f"mismatched {key} in input headers: {known}")
    recorded = [float(c["metadata"]["alphas(rts)"]) for c in curves.values()
                if "alphas(rts)" in c["metadata"]]
    if any(float(c["metadata"].get("p", "-1")) != -1 for c in curves.values()):
        raise ValueError("these observable-correction plots require -p -1")
    if alphas is None:
        if not recorded:
            raise ValueError("no alphas(rts) header; specify --alphas")
        alphas = recorded[0]
    if any(not math.isclose(alphas, value, rel_tol=1e-8, abs_tol=1e-12)
           for value in recorded):
        raise ValueError("--alphas differs from input alphas(rts)")
    if not math.isfinite(alphas) or alphas <= 0:
        raise ValueError("alphas must be finite and positive")
    return alphas


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("hard", "single", "double", "obs", "ll"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--alphas", "--as", type=float,
                        help="reference alpha_s; default: input alphas(rts) header")
    parser.add_argument("--out", type=Path, default=Path("corrections_lambda.png"))
    parser.add_argument("--csv", type=Path, help="default: same stem as --out")
    parser.add_argument("--no-plot", action="store_true", help="write CSV only")
    args = parser.parse_args(argv)
    try:
        curves = {name: read_curve(getattr(args, name), correction=name == "obs")
                  for name in ("hard", "single", "double", "obs", "ll")}
        alphas = check_parameters(curves, args.alphas)
        ratios = build_ratios(curves, alphas)
        csv_path = args.csv or args.out.with_suffix(".csv")
        input_paths = {getattr(args, n).resolve() for n in curves}
        if csv_path.resolve() in input_paths or (
                not args.no_plot and args.out.resolve() in input_paths):
            raise ValueError("output path would overwrite an input file")
        if not args.no_plot and csv_path.resolve() == args.out.resolve():
            raise ValueError("--out and --csv must be different paths")
        if not args.no_plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(ratios)
            writer.writerows(zip(*ratios.values()))
        print(f"Saved {csv_path} (alpha_s = {alphas:g})")
        skipped = sum(ll <= 0 for ll in curves["ll"]["y"])
        if skipped:
            print(f"Warning: {skipped} nonpositive-LL points masked as NaN.",
                  file=sys.stderr)
        if not curves["ll"]["metadata"]["_sl_observable"] and any(
                curves[n]["metadata"]["_sl_observable"]
                for n in ("hard", "single", "double")):
            print("Warning: physical LL versus SL components; preserving the "
                  "requested ratios, which also include a measurement difference.",
                  file=sys.stderr)
        if not args.no_plot:
            fig, ax = plt.subplots(figsize=(7.2, 4.8))
            for name, label in (
                    ("hard", r"$(\mathrm{hard}-\mathrm{LL})/(\alpha_s\,\mathrm{LL})$"),
                    ("single", r"$(\mathrm{single}-\mathrm{LL})/(\alpha_s\,\mathrm{LL})$"),
                    ("double", r"$(\mathrm{double}-\mathrm{LL})/(\alpha_s\,\mathrm{LL})$"),
                    ("obs", r"$\mathrm{obs}/(\alpha_s\,\mathrm{LL})$")):
                ax.plot(ratios["lambda"], ratios[name], "o-", label=label,
                        markersize=4, linewidth=1.5)
            ax.axhline(0, color="0.5", linewidth=0.7)
            ax.set(xlabel=r"$\lambda=\alpha_s\ln(Q/E_T)$",
                   ylabel=r"Correction / $(\alpha_s\,\mathrm{LL})$")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=9)
            fig.tight_layout()
            args.out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(args.out, dpi=180)
            plt.close(fig)
            print(f"Saved {args.out}")
    except ImportError:
        parser.error("plotting requires Matplotlib: python3 -m pip install "
                     "matplotlib (or use --no-plot)")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())

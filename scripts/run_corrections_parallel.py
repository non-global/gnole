#!/usr/bin/env python3
"""Launch the six GNole correction samples concurrently.

Defaults reproduce the requested common settings:
  -nev 1e7 -drap 2.0 -nbins 40 -xmur 1.0 -xQ 1.0 -as 0.1
  -lnktmax 10 -etamax 3 -seed 1 -mc-errors

Each process is single-threaded. The OS schedules them across available cores;
this is not fixed CPU affinity. All six start before this script waits.

Modes:
  hard:           -order 1 -hard-only -sl-obs
  single:         -order 1 -single-insertion-only -sl-obs
  double_gauge:   -order 1 -double-insertion-gauge-only -sl-obs
  double_fermion: -order 1 -double-insertion-fermion-only -sl-obs
  obs:            -order 0 -obs-correction-only
  ll:             -order 0

--ll-sl adds -sl-obs only to LL, if a common SL baseline is wanted.
Each job writes NAME.dat and NAME.log. commands.json records the exact commands
and exit codes. Existing targets are protected unless --overwrite is supplied.
Interrupting the launcher stops its children; partial checkpoints/logs are kept.

Usage:
  python3 run_corrections_parallel.py --dry-run
  python3 run_corrections_parallel.py --out-dir correction_runs
  python3 run_corrections_parallel.py --nev 1000 --out-dir smoke_runs

This launcher needs only the Python standard library. It does not rebuild GNole
or plot automatically; after success it prints a plot_corrections.py command.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import time


MODES = {
    "hard": ["-order", "1", "-hard-only", "-sl-obs"],
    "single": ["-order", "1", "-single-insertion-only", "-sl-obs"],
    "double_gauge": ["-order", "1", "-double-insertion-gauge-only", "-sl-obs"],
    "double_fermion": ["-order", "1", "-double-insertion-fermion-only", "-sl-obs"],
    "obs": ["-order", "0", "-obs-correction-only"],
    "ll": ["-order", "0"],
}


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def events(value):
    number = positive(value)
    if number != math.floor(number) or number > 2147483647:
        raise argparse.ArgumentTypeError("must be an integer from 1 to 2147483647")
    return value  # Keep scientific notation in the recorded command.


def stop_children(processes):
    for process in processes.values():
        if process.poll() is None:
            process.terminate()
    for process in processes.values():
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


def main(argv=None):
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gnole", type=Path, default=script_dir.parent/"bin/gnole")
    parser.add_argument("--out-dir", type=Path, default=Path("correction_runs"))
    parser.add_argument("--nev", type=events, default="1e7")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--alphas", "--as", type=positive, default=0.1)
    parser.add_argument("--drap", type=positive, default=2.0)
    parser.add_argument("--nbins", type=int, default=40)
    parser.add_argument("--xmur", type=positive, default=1.0)
    parser.add_argument("--xq", "--xQ", type=positive, default=1.0)
    parser.add_argument("--lnktmax", type=positive, default=10.0)
    parser.add_argument("--etamax", type=positive, default=3.0)
    parser.add_argument("--ll-sl", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the six commands without writing or launching")
    parser.add_argument("--overwrite", action="store_true",
                        help="allow replacement of this launcher's output/log files")
    args = parser.parse_args(argv)
    if args.nbins < 1:
        parser.error("--nbins must be positive")
    binary, outdir = args.gnole.resolve(), args.out_dir.resolve()
    common = [
        "-nev", args.nev, "-drap", str(args.drap), "-nbins", str(args.nbins),
        "-xmur", str(args.xmur), "-xQ", str(args.xq), "-as", str(args.alphas),
        "-lnktmax", str(args.lnktmax), "-etamax", str(args.etamax),
        "-seed", str(args.seed), "-mc-errors",
    ]
    commands = {}
    for name, mode in MODES.items():
        extra = ["-sl-obs"] if name == "ll" and args.ll_sl else []
        commands[name] = [str(binary)] + common + mode + extra + [
            "-out", str(outdir/(name+".dat"))]
    for command in commands.values():
        print(shlex.join(command), flush=True)
    if args.dry_run:
        return 0
    targets = [outdir/(n+suffix) for n in MODES for suffix in (".dat", ".log")]
    manifest = outdir/"commands.json"
    targets.append(manifest)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        parser.error(f"GNole executable not found or not executable: {binary}")
    if outdir.exists() and not outdir.is_dir():
        parser.error(f"output directory is not a directory: {outdir}")
    if any(p.resolve() == binary for p in targets):
        parser.error("an output path would overwrite the GNole executable")
    existing = [p for p in targets if p.exists() or p.is_symlink()]
    if existing and not args.overwrite:
        parser.error("targets already exist; choose a new --out-dir or use "
                     "--overwrite: " + ", ".join(map(str, existing)))
    if any(p.is_dir() for p in existing):
        parser.error("an output target is a directory")

    record = {"commands": commands, "ll_measurement": "SL" if args.ll_sl else "physical",
              "returncodes": {}, "interrupted": False}
    processes = {}
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

    def interrupted(signum, frame):
        raise KeyboardInterrupt

    previous_term_handler = signal.signal(signal.SIGTERM, interrupted)
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps(record, indent=2)+"\n")
        # Launch every child before entering the wait loop.
        for name, command in commands.items():
            with (outdir/(name+".log")).open("w") as log:
                processes[name] = subprocess.Popen(
                    command, stdout=log, stderr=subprocess.STDOUT, env=env)
            print(f"Started {name}: PID {processes[name].pid}", flush=True)
        pending = set(processes)
        while pending:
            for name in list(pending):
                code = processes[name].poll()
                if code is not None:
                    pending.remove(name)
                    record["returncodes"][name] = code
                    print(f"Finished {name}: exit {code}", flush=True)
            if pending:
                time.sleep(0.25)
    except (KeyboardInterrupt, OSError) as error:
        stop_children(processes)
        record["interrupted"] = isinstance(error, KeyboardInterrupt)
        record["returncodes"] = {n: p.returncode for n, p in processes.items()}
        print("Stopped; partial outputs and logs retained." if record["interrupted"]
              else f"Launch failed: {error}", file=sys.stderr)
        return 130 if record["interrupted"] else 1
    finally:
        signal.signal(signal.SIGTERM, previous_term_handler)
        if outdir.is_dir():
            manifest.write_text(json.dumps(record, indent=2)+"\n")
    if any(code != 0 for code in record["returncodes"].values()):
        print(f"One or more jobs failed; check the logs in {outdir}", file=sys.stderr)
        return 1
    plot_command = [sys.executable, str(script_dir/"plot_corrections.py")]
    for name in MODES:
        plot_command += ["--"+name.replace("_", "-"), str(outdir/(name+".dat"))]
    plot_command += ["--alphas", str(args.alphas), "--out",
                     str(outdir/"corrections_lambda.png")]
    print("All six runs completed. Plot with (requires Matplotlib):\n"+
          shlex.join(plot_command))
    return 0


if __name__ == "__main__":
    sys.exit(main())

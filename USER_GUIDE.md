# GNole user guide

This guide provides some guidance with the use of gnole-v2.0.
Updated in October 2026.

The supported accuracies are LL (`-order 0`) and expanded NLL (`-order 1`).
The unexpanded NLL implementation and the old `-expand-nll` option are removed.
Omitting `-order` selects LL, not NLL.

## Quick start

### LL with the physical measurement

```sh
./bin/gnole -order 0 -nev 1e6 -seed 1 \
  -drap 2 -nbins 40 -as 0.118 -xmur 1 -xQ 1 \
  -lnktmax 10 -etamax 5 -mc-errors -out ll.dat
```

### Full expanded NLL with the physical measurement

```sh
./bin/gnole -order 1 -nev 1e6 -seed 1 \
  -drap 2 -nbins 40 -as 0.118 -xmur 1 -xQ 1 \
  -lnktmax 10 -etamax 5 -mc-errors -out nll.dat
```

This includes the hard contribution and the single- and double-insertion
evolution pieces.

### NLL with adaptive insertion sampling

```sh
./bin/gnole -order 1 -nev 1e6 -seed 1 \
  -drap 2 -nbins 40 -lnktmax 10 -etamax 5 \
  -mc-errors -adaptive-insertion -out nll_adaptive.dat
```

Adaptive sampling may improve distribution-bin errors while increasing the
uncertainty of the integrated weight. Compare the errors in the quantity you
actually need; it is not a universal efficiency improvement.

### NLO matching

```sh
./bin/gnole -order 1 -matching -nev 1e6 -seed 1 \
  -drap 2 -nbins 40 -lnktmax 10 -mc-errors -out nll_matching.dat
```

### Soft observable correction only

```sh
./bin/gnole -order 0 -obs-correction-only -nev 1e6 -seed 1 \
  -drap 2 -nbins 40 -lnktmax 10 -etamax 5 \
  -mc-errors -out obs_correction.dat
```

The flag automatically enables the SL measurement. The two output tables
already contain the cumulant correction, with coupling and bin-width factors
included. Do not cumulatively sum them or multiply by another coupling.

## Command-line syntax

GNole uses single-dash, case-sensitive names: for example `-xmur` and `-xQ`.
Supply a value after value-taking options, as in `-nev 1e6`; use a space rather
than an equals sign. Flags are enabled by their presence and take no value.
Omit a flag to disable it: `-sl-obs 0` does not disable SL measurement.

`-nev` is read as a floating-point value but must be a whole event count, so
`-nev 1e6` is valid. Integer-typed options such as `-seed`, `-nbins` and `-order`
require integer tokens: write `-nbins 100`, not `-nbins 1e2`.

If a value-taking option occurs repeatedly, its last value wins. Argument-file
contents are appended and can therefore override command-line values. Unknown
option names and malformed required numbers are rejected in normal runs.

## Complete option reference

There are 28 application options, plus the help function and argument-file
function supplied by the parser. Help accepts both `-h` and `--help`.
All ten application flags default to off.

### Events and output

| Option | Default | Meaning and accepted values |
|---|---|---|
| `-nev N` | `1e6` | Number of complete outer Monte Carlo events. Must be a finite positive whole number fitting the program's `int` counter; currently at most 2,147,483,647 on the supplied platform. Scientific notation is accepted. Use at least two events for defined statistical errors. |
| `-seed S` | `0` | Integer seed for the fixed GSL MT19937 generator. Zero is deterministic, not a request for a time-based seed. Negative integer seeds are accepted and converted to the generator's unsigned seed type; use nonnegative identifiers for clear, portable studies. |
| `-out FILE` | `output.dat` | Output file, overwritten at checkpoints and at completion. The parent directory must already exist and be writable. Quote paths containing spaces in the shell. `-out ''` sends the final output to standard output; `-out -` instead names a literal file `-`. |
| `-mc-errors` | off | Append complete-event Monte Carlo standard errors. Signed partners within an outer event are combined before squaring. These are statistical errors, not perturbative scale bands. |
| `-output-precision D` | `17` | Significant digits for histogram and diagnostic output. Integer from 1 through 17. The parameter comments assembled earlier can still use ordinary stream precision; the recorded command line preserves the requested tokens. |
| `-diagnostics` | off | Add comment lines with emitter, antenna-cache, zero-work, rapidity-repair and RNG counters. Does not add numeric histogram columns. |

### Accuracy and measurement

| Option | Default | Meaning and accepted values |
|---|---|---|
| `-order K` | `0` | `0`: LL. `1`: expanded NLL. No other order is supported. |
| `-drap Y` | `1.0` | Positive full rapidity width of the measurement slice. Membership is strict: `abs(rapidity relative to event axis) < Y/2`. This is not its half-width. |
| `-as A` | `0.118` | Positive strong coupling at the hard scale `Q = sqrt(s)`. It is not an automatically evolved value at another reference scale. Energies in the code use the normalized hard scale; there is no separate collision-energy CLI option. |
| `-xmur R` | `1.0` | Positive renormalization-scale factor, `mu_R = xmur * Q`. The effective coupling computed at this scale must remain finite and positive. |
| `-xQ QF` | `1.0` | Positive resummation-scale factor. The sampling domain requires `log(xQ) < 1/(2*alpha_s(xmur)*b0)`. It is independent of `-xmur`. |
| `-lnktmax LMAX` | `20.0` | Positive infrared cutoff on the evolution variable `ln_kt(t)`. Fractional values are retained. Evolution also has Landau/evolution-time and, at order one, legacy endpoint limits; this option is not the sole stopping condition. |
| `-etamax ETA` | `5.0` | Positive technical rapidity buffer used in the emitting dipole frame. It is not the physical slice width or a fixed laboratory pseudorapidity cut. Check convergence when changing it. |
| `-sl-obs` | off | Use the SL measurement `v_SL = exp(-ell)`, with the already resummation-scale-adjusted log supplied by the shower, instead of physical transverse energy. Can be used with LL, NLL and matching. Automatically enabled by `-obs-correction-only`. |
| `-p P` | `-1.0` | `-1` gives the ordinary plotted log `L = ln(1/v)`. Any positive value gives the modified plotted log `L_p = ln(1+v^(-p))/p`. Zero and other negative values are rejected. Observable-correction-only currently requires `-1`. |

Use finite, physically sensible scale, coupling, cutoff and slice parameters.
Acceptance by the parser is not a validation of an exotic physics setup.

### NLL component selection and matching

| Option | Default | Meaning and compatibility |
|---|---|---|
| `-hard-only` | off | Select LL evolution with the one-loop two-parton hard factor and the three-parton hard contribution; omit single/double NLL evolution insertions. Requires `-order 1`. Accepts `-matching` and `-sl-obs`. |
| `-single-insertion-only` | off | Select the two-parton single-insertion evolution contribution, retaining its LL baseline. Requires `-order 1`; cannot be combined with matching or another component selector. |
| `-double-insertion-only` | off | Select double-insertion gauge and fermion contributions and their counterterms, retaining the LL baseline through the explicit branch-zero evolution. Requires `-order 1`; cannot be combined with matching or another component selector. |
| `-matching` | off | Enable the hard-sector matching prescription. Requires `-order 1`, either full NLL or `-hard-only`. Both physical and SL measurements are accepted. Not permitted with single-only, double-only or observable-correction-only. |
| `-obs-correction-only` | off | Run unchanged LL trajectories and output only the pure NLL soft measurement correction `DeltaSigma = Sigma_full - Sigma_SL`. Requires `-order 0` and `-p -1`; implies SL measurement. Rejects matching, all NLL component selectors and adaptive insertion. |

Use at most one of the three NLL component selectors. With none selected,
`-order 1` runs the full expanded NLL calculation.

### Sampling and numerical controls

| Option | Default | Meaning and scope |
|---|---|---|
| `-variance-mode MODE` | `coupled` | `legacy`: original separate hard counter-events and independent continuations. `merge`: combine identical two-parton hard counterterms. `coupled`: also correlate real/subtraction continuations. These change sampling, not the intended estimator. They affect NLL branches, not the ordinary LL shower. |
| `-adaptive-insertion` | off | Learn a first-insertion importance proposal from completed signed two-parton events. Effective for full NLL and single/double insertion runs. Accepted but performs no insertion learning in ordinary LL or hard-only runs. Rejected in observable-correction-only mode. |
| `-insertion-power N` | `5.0` | Finite first-insertion importance-map power greater than 2.5. Its Jacobians preserve the canonical estimator. Does not affect LL or hard-only runs. |
| `-second-power N` | `1.0` | Finite power greater than 0.5 for the weighted second-insertion map. Relevant only when double-insertion branches are executed; has no effect on LL, hard-only or single-only runs. |
| `-rng-replay MODE` | `tape` | `tape` or `snapshots`. Two implementations of the same correlated uniform replay with the longest-branch endpoint committed. Used in coupled mode; accepted but has no replay role in legacy, merge or ordinary LL evolution. |
| `-reference-kernels` | off | Debugging comparison mode: restore preceding arithmetic and disable the newer speed shortcuts. It does not automatically select legacy variance sampling or six-digit printing. Not a recommended production-precision setting. |


### Histogram controls

| Option | Default | Meaning and scope |
|---|---|---|
| `-nbins N` | `100` | Positive integer number of ordinary L and evolution-time bins. Does not change the independent lambda grid. Correction-only L output has `N+1` shifted bins. |
| `-obsmax X` | value of `-lnktmax` | Positive upper display limit of the ordinary L histogram. This default is dynamic: `-lnktmax 10` without `-obsmax` gives a display maximum of 10. Does not change the shower cutoff. |

The grids are:

| Quantity | Ordinary grid | Observable-correction-only grid |
|---|---|---|
| L | `[0, obsmax]`, `nbins` bins | `[-deltaL/2, obsmax+deltaL/2]`, `nbins+1` bins |
| lambda | `[0, 0.5]`, 10 bins | `[-0.025, 0.525]`, 11 bins |
| Evolution time t | `[0, 1]`, `nbins` bins | Not printed |

Here `deltaL = obsmax/nbins`, and lambda is `as * plotted_L`.
Changing `-as` changes the mapping between L and lambda, but not the fixed
lambda bin boundaries. The time-grid upper bound stays at the compiled
`EVOLCUT=1`; it is not reset to the time corresponding to `-lnktmax`.
Values below the lower edge enter underflow; values at or above the upper
edge enter overflow.

### Help and argument files

| Option | Default | Meaning |
|---|---|---|
| `-h` or `--help` | off | Print built-in help and exit without generating events. Use the flag alone for an overview; parameter validation can occur before help is printed if other arguments are supplied. |
| `-argfile FILE` | none | Read additional whitespace-separated tokens from a file. Supported by the parser even though absent from the ordinary `-h` listing. See the precedence rules below. |

## Argument files

For example, save the following as `nll.args`:

```text
// Expanded NLL at central scales
-order 1
-nev 1e6
-seed 1
-drap 2
-nbins 40
-as 0.118
-xmur 1
-xQ 1
-lnktmax 10
-etamax 5
-mc-errors
-out nll.dat
```

Run it with:

```sh
./bin/gnole -argfile nll.args
```

The file contents are appended after the command-line tokens. Therefore a
file's `-seed 1` overrides an earlier CLI `-seed 2`, even if `-argfile` appeared
before `-seed 2` in the shell command. Avoid repeated options across the two
sources if that precedence is not intended.

Argument files are not shell scripts: there is no shell variable expansion,
quote removal, escaped-space handling or need for backslash continuations.
Tokens containing `//` start a comment that discards the rest of their line;
use `//` comments, not `#` comments. Quoted multiword filenames do not work
as a single token inside these files. Keep argument files simple and avoid
circular file inclusions. Relative filenames resolve from the process's
working directory, not from the argument file's directory.

## Compatibility rules

| Calculation | Order | Component selector | Matching | SL measurement | Adaptive insertion |
|---|---|---|---|---|---|
| Ordinary LL | `0` | None | Rejected | Optional | Accepted, no insertion learning |
| Full expanded NLL | `1` | None | Optional | Optional | Optional |
| Hard sector | `1` | `-hard-only` | Optional | Optional | Accepted, no insertion learning |
| Single insertion | `1` | `-single-insertion-only` | Rejected | Optional | Optional |
| Double insertion | `1` | `-double-insertion-only` | Rejected | Optional | Optional |
| Soft observable correction | `0` | `-obs-correction-only` | Rejected | Automatic | Rejected |

The correction mode additionally requires `-p -1`. Ordinary modes allow the
positive modified-log parameter. Flags such as `-mc-errors`, `-diagnostics`
and `-output-precision` can be used in all these modes.

### Component examples

Use a different output file for each selection:

```sh
./bin/gnole -order 1 -hard-only -nev 1e6 -drap 2 -seed 1 -mc-errors -out hard.dat
./bin/gnole -order 1 -single-insertion-only -nev 1e6 -drap 2 -seed 1 -mc-errors -out single.dat
./bin/gnole -order 1 -double-insertion-only -nev 1e6 -drap 2 -seed 1 -mc-errors -out double.dat
```

For example, `-hard-only -sl-obs -matching` is accepted. Combining hard-only
with single-only or double-only is not.

## Reading ordinary output

Ordinary runs print six blocks:

```text
# dSdt (differential)
# dSdt (cumulative)
# dSdlnET (differential)
# dSdlnET (cumulative)
# dSdlambda (differential)
# dSdlambda (cumulative)
```

The L coordinate uses the chosen ordinary or modified log. In physical mode
the measured value is transverse energy relative to the event axis; in SL
mode it is the supplied evolution-scale measurement. LL evolution fills the
first in-slice emission and stops. NLL branches also make signed insertion
and reconstructed-parent deposits, and physical matching can measure a sum
over the hard partons. This is not an exact sum over a subsequently fully
continued soft shower.

### Columns and normalization

| Block | Without errors | With `-mc-errors` |
|---|---|---|
| Differential | `lower midpoint upper density` | Add column 5: `standard_error` |
| Cumulative | `edge cumulative` | Add column 3: `standard_error` |

For an ordinary bin of width `delta`, the density is its signed weight sum
divided by `nev * delta`. Differential rows ascend in coordinate; cumulative
rows descend. For plotting densities, use midpoint column 2 and density
column 4, with error column 5 when present. For cumulants, use columns 1 and 2,
with error column 3.

The cumulative blocks are inverse sums of recorded histogram fills: each
edge includes overflow and all physical bins at or above that edge.
Underflow is excluded from those numeric rows and added only in the final
`# with_underflow` comment. That final value is total accepted fill weight
per outer event, not a quantity required to equal one.

Histories with no measured hit before the cutoffs are not entered as an
explicit zero-observable atom. They still count in `nev` and in the error
denominator. Therefore the printed recorded-fill cumulant should not be
identified automatically with a unit-normalized probability including the
unresolved tail. NLL and component weights can be signed; do not clip negative
entries merely to impose probability-like behavior.

### Flow bins

Ordinary differential flow lines report **raw weight sums**, not divided by
`nev` or by a bin width. Their optional `standard_error_of_mean` is already
the standard error of the per-event flow weight. To obtain a flow-bin mean,
divide its raw central weight by the event count, but do not divide that
already normalized error by `nev` again.

The ordinary overflow line starts with `over flow bin` without `#`.
Consequently the entire file is not a plain rectangular numeric table:
select the desired block's numeric rows and handle comments and flow lines
separately when loading it into plotting software.

### Statistical errors

All signed contributions from one complete outer event are summed before
its variance contribution is formed. Cumulative errors include covariance
between bins and flow contributions; they are not obtained by independent-bin
quadrature. Zero-fill events are included. At least two complete events are
needed: `-nev 1 -mc-errors` is accepted but gives undefined, NaN errors.

These errors describe Monte Carlo uncertainty of the weighted estimator.
They are not scale-variation bands, and rare signed tails can make finite
error estimates unstable. Use multiple seeds and adequate statistics when
quantifying differences.

## Reading observable-correction output

This mode prints only:

```text
# observable correction (L)
# observable correction (lambda)
```

Each numeric row is:

```text
lower center upper delta_cumulant [standard_error]
```

For the first in-slice LL emission, the sampled log is `ell`, and

\[
c=\log\frac{v_{\rm SL}}{v_{\rm full}}
 =-\ell-\log k_{T,\rm Born}.
\]

The histogram receives `w*c` at `L_SL=ell` and `lambda_SL=as*ell`. Defining
`W_j` as the sum of those signed weights in a correction bin, the printed
values are

\[
\Delta\Sigma_j^{(L)}=\frac{W_j}{N_{\rm ev}\,\Delta L},\qquad
\Delta\Sigma_j^{(\lambda)}=
 \frac{\alpha_s(Q)\,W_j}{N_{\rm ev}\,\Delta\lambda}.
\]

The errors use the same factors. Here the coupling is the input `-as`, not
the running coupling at each emission. All outer events, including no-hit
and zero-correction events, contribute to the normalization. An in-slice
emission with zero correction still terminates the LL shower.

The result is the linearized pure NLL cumulant correction
`Sigma_full - Sigma_SL`, not `DeltaSigma/as` and not an ordinary density to
integrate. No extra coupling, color factor or `1/(2*pi)` should be applied.
It is not the exact finite-coupling difference of two complete distributions.

The shifted bin centers align with ordinary cumulant edges. They are finite
bin averages of the correction function; support-boundary rows are not exact
endpoint values. L and lambda use independent grids, so their rows need not
agree pointwise without corresponding physical bin boundaries. Refine the L
grid with `-nbins` and check cutoff dependence.

Correction flow comments, unlike ordinary differential flow comments, already
report raw weighted sums **per event**, `W_flow/nev`, with unscaled per-event
errors. They are labelled `not DeltaSigma` because no correction density
normalization is assigned to a flow bin.

Do not add this correction to an already physical-measurement result without
checking for double counting. It also does not restore physical hard-event
matching to SL endpoint deposits. If combining LL and correction estimates
from shared random histories, include their covariance; the current output
does not supply cross-run covariance. Independent-run errors can be combined
only when the independence assumption is appropriate.

## Practical running advice

### Scale variations

Use separate filenames and hold the other inputs fixed:

```sh
./bin/gnole -order 1 -nev 1e6 -drap 2 -nbins 40 -lnktmax 10 -seed 1 \
  -xmur 0.5 -xQ 1 -mc-errors -out nll_mur05_q1.dat
./bin/gnole -order 1 -nev 1e6 -drap 2 -nbins 40 -lnktmax 10 -seed 1 \
  -xmur 1 -xQ 2 -mc-errors -out nll_mur1_q2.dat
```

The program does not construct an envelope or a perturbative uncertainty band
automatically. Matching seeds can help comparisons but do not guarantee identical
histories when parameters or sampling modes change.

### Cutoffs and display range

To compare infrared cutoffs on the same display grid, fix `-obsmax` explicitly,
for example compare `-lnktmax 10 -obsmax 10` with `-lnktmax 12 -obsmax 10`.
Similarly, hold `-drap` fixed while checking the technical `-etamax` buffer.
An unchanged output after increasing one cutoff can mean another stopping
limit is already dominant; it is not, by itself, a complete convergence test.

### Reproducing earlier sampling

For comparisons with the preceding variance-reduced edition, use the same
physics, seed, variance and map settings, plus

```text
-reference-kernels -output-precision 6
```

To additionally use the earlier expanded-only sampling convention, select
`-variance-mode legacy`. For explicit checkpoint replay choose
`-rng-replay snapshots`. These controls do not restore the removed unexpanded
algorithm and do not promise agreement with every bug in the original archive.
Reproducibility also depends on compiler, optimization and library versions.

### Diagnostics and progress

`-diagnostics` adds emitter-call counts, mean and maximum dipole multiplicity
at selection, zero-work skips, antenna evaluations/cache hits, nonfinite-weight
rejections, rapidity-sum rebuilds, actual generator uniform draws, snapshot
bytes copied and tape capacity. Tape capacity is counted in stored doubles,
not bytes. Tape can generate fewer actual uniforms than snapshots while
committing the same future stream. Zero-work skips are branch/continuation
counts, not a count of empty outer events.

Adaptive runs also print proposal update counts and probabilities. Progress
messages go to standard error. Checkpoints initially occur every million
events, with the interval increasing for long runs; there is no CLI option
for that interval. Checkpoints overwrite the same output file, and a final
write is always performed. They are not resumable generator-state snapshots.

Check the output directory before a long run: the output file is opened at a
checkpoint or final write, not necessarily when the calculation starts.

### Errors and removed options

Invalid values and incompatible combinations exit with status 1; help exits
with status 0. Input validation occurs before event generation and output
writing, while runtime failures can leave an earlier checkpoint. Typical
examples are mutually exclusive selectors, matching at order zero,
observable correction at order one, and missing values for numeric options.

Do not pass `-expand-nll`; it is unrecognized in this edition. There are no
CLI options for changing the color factors, flavor count, lambda binning,
evolution-time limit, replay generator type or output interval. These are
compiled settings, not undocumented switches.

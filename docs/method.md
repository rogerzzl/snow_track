# Method: multi-layer, dual-phase source tracking

This page states the scheme implemented in `src/snowtrack/tracker.py` and maps
every equation to the code. It follows Supplementary Text S3 of Zhu et al.
("Subsurface-routed snowmelt buffers river response to warming on the Tibetan
Plateau"); the generalisations (any number of layers, any set of sources,
optional ice phase, unexplained water) are marked as such.

## 1. What is tracked

A model cell has a soil column of $N \ge 1$ layers (top to bottom
$n = 1 \dots N$). The host model supplies, per time step, the bulk liquid
storage $W_n^{liq}(t)$ (and optionally the ice storage $W_n^{ice}(t)$) of every
layer, the bulk fluxes between layers and to the channel, and the water that
reaches the ground from each source $k$ (rain, seasonal snowmelt, glacier
melt, ... any set). The tracker carries source-resolved storages
$W^{liq}_{n,k}$, $W^{ice}_{n,k}$ for every layer and every source $k \in
\{\text{sources}\} \cup \{U\}$, where $U$ ("UNKNOWN") holds water whose origin
is not resolved (initial storage, and fluxes the model did not export), such
that at every step

$$\sum_k W^{liq}_{n,k}(t) = W^{liq}_n(t), \qquad \sum_k W^{ice}_{n,k}(t) = W^{ice}_n(t).$$

Initialisation (`SourceTracker.run`, `t == 0`): $W_{n,U}(t_0) = W_n(t_0)$ for
both phases, all other sources zero; runoff of the first step is assigned to
$U$. The scheme is valid without spin-up: $U$ simply reports the share of
runoff originating from antecedent water, and it decays as the stores turn
over (`FW_UNKNOWN`, `RUNOFF_UNKNOWN_*` in the output).

## 2. Surface boundary

$$S_{tot} = \sum_k S_k, \qquad f_k = S_k / S_{tot} \ (S_{tot} > \epsilon),$$
$$i_k = I\, f_k, \qquad R_{surf,k} = R_{surf}\, f_k,$$

i.e. surface water is well mixed on the infiltration time scale and
infiltration does not select a source (code: `frac`, `inflow`, `surf[t]`).

*Generalisation – unexplained water.* Many models export direct runoff that on
some days exceeds the supply that reached the ground (e.g. runoff generated
from storage, lake or wetland terms), or infiltration not matched by the
supply. The tracker attributes only the part that the supply can explain,
$\min(R_{surf}, S_{tot} - \min(I, S_{tot}))$ and $\min(I, S_{tot})$, by the
source fractions; the remainder goes to $U$ so that
$\sum_k R_{surf,k} = R_{surf}$ and $\sum_k i_k = I$ always hold
(`RUNOFF_UNKNOWN_SURFACE`). With `--strict-surface-closure` a violation of
$S_{tot} = I + R_{surf}$ raises instead. For the VIC output of the study this
rule reproduces the archived surface attribution exactly.

## 3. Transport and evapotranspiration (proportional partitioning)

For layer $n$ the model supplies percolation $Q^{perc}_n$ (to layer $n+1$),
drainage/baseflow to the channel $Q^{base}_n$, and evapotranspiration
withdrawn from the layer $E^{tr}_n$. Complete mixing within the layer's liquid
store defines the composition

$$c^{liq}_{n,k} = W^{liq}_{n,k} / W^{liq}_n \quad (W^{liq}_n > \epsilon),$$

and every liquid outflux carries it:

$$q^{perc}_{n,k} = Q^{perc}_n c^{liq}_{n,k}, \quad q^{base}_{n,k} = Q^{base}_n c^{liq}_{n,k}, \quad e_{n,k} = E^{tr}_n c^{liq}_{n,k}, \quad q^{in}_{n,k} = \begin{cases} i_k & n = 1 \\ q^{perc}_{n-1,k} & n > 1.\end{cases}$$

Summing over $k$ recovers the bulk fluxes by construction. Percolation out of
the bottom layer is a loss from the column (deep drainage). If a layer is
empty ($W^{liq}_n < \epsilon$) but the model still withdraws water from it,
the withdrawal is attributed to $U$ (closure is kept; the budget residual
diagnostics flag the inconsistency).

**Which composition carries the withdrawals** (`TrackerConfig.composition`):

| option | withdrawals use | reported baseflow uses | note |
|---|---|---|---|
| `end` (default) | composition of the store *after* the step's inflow has been mixed in, $c = (W + q^{in}) / \sum(W + q^{in})$ | same | implicit well-mixed reservoir; inflow of a step can leave in the same step |
| `start` | composition at the start of the step | same | explicit; version 1 of this repository |
| `legacy` | composition at the start of the step | composition at the end of the step, after reconciliation (section 5) | convention of the code that produced the archived results of Zhu et al. |

All three conserve mass exactly; they differ only in the timing of mixing
within one step. For the daily VIC-glacier cells of the study the annual
subsurface snowmelt differs by less than 2 % between conventions
(`tests/test_regression_vic.py`).

## 4. Phase change (dual-phase coupling)

With bulk freezing $\Phi^{freeze}_n$ (liquid to ice) and melting
$\Phi^{melt}_n$ (ice to liquid), and the ice composition
$c^{ice}_{n,k} = W^{ice}_{n,k} / W^{ice}_n$:

$$\phi^{freeze}_{n,k} = \Phi^{freeze}_n c^{liq}_{n,k}, \qquad \phi^{melt}_{n,k} = \Phi^{melt}_n c^{ice}_{n,k}.$$

Phase change is non-fractionating with respect to source, so
$S_{n,k} = W^{liq}_{n,k} + W^{ice}_{n,k}$ is unchanged by it: soil ice stores
the historical source mixture and releases it later ("memory"). If
$W^{liq}_n < \epsilon$ no freezing is applied, if $W^{ice}_n < \epsilon$ no
melting (near-zero guards, S3 step 14). When a model exports ice storage but
no phase-change fluxes, the net phase change is inferred from the change of
the ice storage (`CellData.from_dataframe(infer_phase_change=True)`).

## 5. Governing equations, operator splitting and reconciliation

$$\frac{dW^{liq}_{n,k}}{dt} = q^{in}_{n,k} - q^{perc}_{n,k} - q^{base}_{n,k} - e_{n,k} + \phi^{melt}_{n,k} - \phi^{freeze}_{n,k}, \qquad \frac{dW^{ice}_{n,k}}{dt} = \phi^{freeze}_{n,k} - \phi^{melt}_{n,k}.$$

Each step applies (A) transport + ET with the proportional partitioning, then
(B) freeze–thaw with the phase compositions (first-order explicit splitting,
S3 step 13). After each sub-step tiny negative values from round-off are set
to zero; larger negative values mean the exported fluxes withdraw more than
the store holds – they are clipped and the clipped mass is reported
(`CLIPPED_MASS`), or an error is raised with `--strict`.

**Reconciliation (S3 step 16).** The host model's bulk storage is
authoritative. At the end of every step the tracked storages of each layer
and phase are rescaled so that their sum equals the exported bulk storage.
The difference before rescaling, `BUDGET_RESIDUAL_LIQ<n>` / `_ICE<n>`, is the
part of the model's own water balance that its exported fluxes do not
explain; it is zero for a consistent export (Example 01, the synthetic cells)
and small but non-zero for models that do not export every term. If the tracked store is empty while
the bulk store is not, the bulk amount is assigned to $U$
(`UNKNOWN_INJECTED`). In the `legacy` convention the reconciliation is the
literal normalisation of the archived code (fractions clipped to $[0, 1]$ and
renormalised), which is why it reproduces the archive bit for bit.

## 6. Surface vs subsurface snowmelt runoff

$$R_{surf,snow}(t) = R_{surf}(t)\, f_{snow}(t), \qquad B_{snow}(t) = \sum_n q^{base}_{n,snow}(t),$$

(`RUNOFF_SNOW_SURFACE`, `RUNOFF_SNOW_SUBSURFACE`; `RUNOFF_SNOW_TOTAL` is their
sum; analogous columns for every other source and for $U$;
`F_<SRC>_<PATH>` are the shares of `TOTAL_RUNOFF` $= R_{surf} + \sum_n Q^{base}_n$).
`RUNOFF_SOURCE_CLOSURE_ERROR` $= \sum_k (R_{surf,k} + B_k) - \text{TOTAL\_RUNOFF}$
is zero to round-off in every configuration.

## 7. Reduction to one store

With $N = 1$ and no ice columns the scheme is a single well-mixed soil store:
$W_1$ = the exported soil-water storage, $Q^{base}_1$ = baseflow,
$E^{tr}_1$ = evapotranspiration withdrawn from the soil,
$I = \max(S_{tot} - R_{surf}, 0)$ when no infiltration term is exported.
This is the minimum any model can supply; every layered or dual-phase term
that a model exports is used by the same code without any change
(Example 01, `configs/examples/clm_netcdf.yaml`). `--composition legacy`
reproduces the archived per-cell files of the study to machine precision
(`docs/reproduce_paper.md`).

## 8. Metrics (Methods and S4, S7, S8, S9)

Implemented in `metrics.py` and `lagdiag.py` on an outlet or cell series:

* **Fractions**: surface / subsurface snowmelt divided by total runoff at the
  same resolution (daily, monthly climatology, annual, overall).
* **T50**: first day of the (melt-)year on which the cumulative snowmelt
  runoff reaches 50 % of the annual total; **centroid** $= \sum t\,q(t) / \sum q$.
  Calendar years for the trends of Table S4; melt-response years (1 April to
  31 March) for the lag diagnostics.
* **Timing trend**: OLS slope of annual T50, reported in days per decade with
  advance positive.
* **Buffering index** $B = \text{trend}_{surf} - \text{trend}_{Q}$ (S7).
* **Late-season subsidy** (S8): July–October subsurface snowmelt fraction of
  seasonal runoff (years with seasonal total < 1 mm excluded) and the minimum
  30-day moving-average runoff of the window; early (1968–1988) vs late
  (1998–2018) means and their differences.
* **Apparent lag** (S9): non-negative distributed-lag model
  $Q_{sub}(t) = \sum_j \beta_j \overline{Q_{surf}}_{[\tau \in B_j]}(t-\tau)$ with bins
  0–30, 30–90, 90–180, 180–365, 365–730, 730–1825 days, uniform daily kernel
  within a bin, apparent mean / median / 90th-percentile lag, long-lag
  fractions, apparent origin-year allocation, annual T50 delay
  $T_{50,sub} - T_{50,surf}$ and the next-winter release ratio. These are
  effective storage-release time scales of the model, not water ages.

## 9. Routing to the basin outlet

`snowtrack route` (`routing.py`) applies the Lohmann channel routing scheme
(Lohmann et al., 1996, 1998): the tracked components of every cell are
convolved with the impulse response of the linearised Saint-Venant equation
(diffusion wave) over the flow distance $x$ from the cell to the gauge,

$$h(x, t) = \frac{x}{2\sqrt{\pi D t^3}} \exp\!\left(-\frac{(x - C t)^2}{4 D t}\right),$$

with wave velocity $C$ and diffusivity $D$. Surface components use the
surface-pathway velocity, subsurface components the subsurface-pathway
velocity; routed totals are the sums of the routed components. $h$ is the
inverse-Gaussian density with mean $x/C$ and shape $x^2/(2D)$, so the mass
arriving in each time step is taken exactly from its cumulative distribution.
**The two velocities and the diffusivity are calibration parameters**,
calibrated against gauged discharge together with the hydrological model; the
command requires them explicitly. Flow distances come from the cell
coordinates (great-circle distance to the gauge) or from a cells table with
flow-path distances or prescribed lags of a river-network routing
(`--cells`). A pure translation scheme (`--scheme lag`) is available for
comparison.

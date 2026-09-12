# Reporting retained tool failures

Added during execution on 2026-09-12, after the first failure was observed.
This is a reporting addendum, not part of the original pre-execution freeze.

Shepherd stopped on random library seed `20260912343`: its automatic
substitution-error estimator could not produce an acceptable estimate.
The tool requests a supplied error-rate estimate. The registered invocation
uses automatic estimation, so this failure is retained with its original
command and logs. No tool parameter, seed, input, or timed invocation is
replaced. Additional failures, if any, receive the same treatment.

The original frozen `analyze.py` requires all cells to succeed. It remains
unchanged. `analyze_available.py` makes the incomplete evidence reviewable:

- Every registered cell remains in the per-input export, including failures.
- Every method/design reports successes and failures. Full-design means and
  medians are unavailable if any of its registered results are missing.
- Successful-only descriptive summaries are exported separately, explicitly
  conditional on successful execution and excluded from leader highlighting.
- A contrast is computed only if **all 60 registered pairs** are present.
  Otherwise it is unavailable: no reduced-sample test or imputation.
- Available contrasts use the original paired function, fixed bootstrap seed,
  thresholds and family size of eight. Unavailable contrasts occupy their
  original positions and contribute a conservative placeholder of one only
  to the internal Holm calculation; they have no reported p-value or claim.
- No missing score is treated as an accuracy loss or win. Successful runtime
  cannot stand in for a failed invocation's time-to-valid-output.

This limits the conclusions we can draw about an incomplete competitor/design
without discarding complete comparisons on the other registered inputs.
Any statistical statement must disclose this reporting addendum and identify
unavailable comparisons. Completion is an additional observed outcome, not
evidence of universal reliability or accuracy superiority.

Background: missing outputs can make comparisons conditional on method
success, and replacement, exclusion or retuning changes their interpretation;
frequency and handling should be reported explicitly. See
[Pawel et al., Handling Missingness, Failures, and Non-Convergence in Simulation Studies](https://arxiv.org/html/2409.18527v3).

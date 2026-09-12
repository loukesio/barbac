# Timing observations affected by system sleep

Recorded during the campaign on 2026-09-12, before any repair invocation.
The Mac power log records idle sleep at 22:44:52 local time and a return to
full wake at 23:08:13 (UTC+02). Several maintenance sleep and partial-wake
intervals occurred between those endpoints. This explains long calendar-time
gaps that do not appear in the recorded monotonic durations.

All original clustering outputs, scores and timings remain preserved.
Accuracy is not replaced. The affected timing block is selected solely by
overlap of receipt timestamps with that power-state interval, including a
one-second margin for the power log's timestamp resolution. Partial-wake
observations are included conservatively. No selection uses a runtime value,
accuracy result, ranking or p-value.

The 28 successful cells in this block will each be measured once more, in
their original relative order, with idle sleep prevented. This is an explicit
exception to the original single-invocation plan, needed because the machine
did not remain awake. It is not an additional accuracy campaign. The five
already failed calls are retained as failures and are not retried.

The same frozen tool code, library, input and settings are used. Each new call
uses the original fresh-worker timing boundary. Its centroid and member files
must match the original files byte-for-byte before its timing is accepted.
Both measurements are exported, with the repair identified; the new valid
observation replaces the affected observation in the displayed time summary.
There is no choosing the faster of the two. An unsuccessful repair remains
unavailable and is not retried automatically.

`timing_repair_protocol.json` records the exact affected cells, original
receipt hashes, source fingerprint and power-log excerpt before retiming.
The repaired records live under `generated/timing_repair/`, separately from
the original campaign. `timing_results.json` and `publication_results.csv`
retain the timing selection explicitly. Accuracy analysis and its
failure-reporting addendum are unchanged.

Original timing summaries remain in `summary.json` and `speed_contrasts.json`.
The publication uses `publication_summary.json` and `timing_contrasts.json`.
Any displayed speed comparison must disclose this repair. All timings are
single observations per accepted input/configuration; the single Milo
reference in particular cannot establish a precise general speed ratio.

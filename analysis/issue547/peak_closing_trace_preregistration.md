# Peak closing-rate trace preregistration

**Question.** What contact and timestep produce the 32.63% increase in peak absolute closing rate from 0.05 to 0.025 Myr, and how close is that step to a topology event?

**Design.** Repeat the matched 50 Myr seed-748935664 replay at 0.05 and 0.025 Myr for baseline commit `f04cbc0772eca6d2484f7477b20004deb368206c` and candidate commit `fb516de9ee61cb54bb29a9a7050a9fa1f6c07481`. Use node density 4, `PYTHONHASHSEED=0`, and run all four cases serially. These replays repeat earlier cases because their saved output retained only the aggregate peak.

**Readout.** Capture the global maximum absolute contact-closing rate, its step time, plate-pair IDs, signed rate, force, capacity, yield status, plate count, and any topology event during that step. From each run's topology-event log, also report the nearest event before and after the peak step, including time and step distance. Compare peak values and locations across resolutions and arms. This is a diagnostic trace, not a new acceptance threshold.

**Interpretation rule.** Treat an event on the peak step or within one timestep as coincident at this resolution; otherwise report the measured distance without attributing causality. A contact ID match across runs does not establish physical identity after topology changes. Report runtime, peak RSS, and ledger closure, and keep full replay outputs in scratch.

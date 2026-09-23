# Fixed-K magma destination search experiment

## Method

On the issue #205 magma model, use seed 0, node and climate density 4.0, fluid density 2.0,
and 100,000-year steps. Capture the pending parcels and eligible continental destination
nodes immediately before transport passes at steps 4, 20, and 40. The eligible-node filter
from PR #224 is applied to both the complete radius search and fixed-K variants. Compare
K = 64, 128, 256, and 512 with the complete 1,000 km radius search on the **same** snapshot
at each step. For each result, run the same per-parcel weight normalization and global
per-node deposit cap as `run_magma_transport`. These are single unprofiled search timings,
not full-step or full-animation timings.

`Node L1` is the sum of absolute per-node deposit differences divided by the complete
search's total deposit. `Region L1` uses the same calculation after summing deposits into
10-degree latitude/longitude cells. Both count mass moved from one bin to another twice;
they measure distribution change, not net lost volume.

## Results

| Step | Parcels | Complete pairs | K | Selected pairs | Search seconds | Placed volume / complete | Node L1 | Region L1 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 23,168 | 5,730,797 | 64 | 959,681 | 0.111 | 0.9971 | 0.569 | 0.288 |
| 4 | 23,168 | 5,730,797 | 128 | 1,844,044 | 0.191 | 0.9995 | 0.298 | 0.159 |
| 4 | 23,168 | 5,730,797 | 256 | 3,428,632 | 0.380 | 1.0000 | 0.080 | 0.046 |
| 4 | 23,168 | 5,730,797 | 512 | 5,257,900 | 0.654 | 1.0000 | 0.006 | 0.004 |
| 4 | 23,168 | 5,730,797 | complete | 5,730,797 | 0.676 | 1.0000 | 0 | 0 |
| 20 | 60,841 | 7,459,732 | 64 | 1,554,578 | 0.165 | 0.9993 | 0.524 | 0.246 |
| 20 | 60,841 | 7,459,732 | 128 | 2,767,119 | 0.292 | 0.9999 | 0.276 | 0.145 |
| 20 | 60,841 | 7,459,732 | 256 | 4,858,767 | 0.560 | 1.0000 | 0.068 | 0.042 |
| 20 | 60,841 | 7,459,732 | 512 | 6,920,611 | 0.966 | 1.0000 | 0.005 | 0.004 |
| 20 | 60,841 | 7,459,732 | complete | 7,459,732 | 1.072 | 1.0000 | 0 | 0 |
| 40 | 88,300 | 7,504,657 | 64 | 1,625,968 | 0.181 | 0.9966 | 0.576 | 0.268 |
| 40 | 88,300 | 7,504,657 | 128 | 2,818,054 | 0.318 | 0.9992 | 0.288 | 0.153 |
| 40 | 88,300 | 7,504,657 | 256 | 4,889,060 | 0.584 | 1.0000 | 0.073 | 0.045 |
| 40 | 88,300 | 7,504,657 | 512 | 7,001,553 | 1.143 | 1.0000 | 0.005 | 0.004 |
| 40 | 88,300 | 7,504,657 | complete | 7,504,657 | 1.213 | 1.0000 | 0 | 0 |

## Decision

Keep the complete radius search as the default. Fixed K is available only by passing the
optional `max_destinations_per_parcel` argument to `run_magma_transport`. K = 64 gives a
large search speedup but substantially changes where melt lands. K = 256 is a plausible
quality/speed compromise for further study, but it still changes roughly 4-5% of the
10-degree regional deposit totals by the L1 measure. The snapshots cover only 4 Myr of
simulation history; this experiment does not establish its long-run land-fraction impact.

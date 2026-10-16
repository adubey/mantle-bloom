# RFC 0011 gate G-R

[RFC 0011](../../rfc/0011-contact-strength.md)'s gate G-R answers two questions that #452
left open:
- **What blocks the overlap:** a craton indenting weaker crust, or craton on both sides of the
  contact?
- **Do jammed pairs come to matching motion on their own?**

Compact numbers are in [gate_gr.json](gate_gr.json): one summary per seed, and `pooled`.

## Method

The runs use seeds 326494039, 748935664, 580014012 and 111065103, at #452's controls:
`generate_world` defaults and 1 Gyr. They run at 2 Myr steps, one at a time, on `6579830`
with this branch's `bin/debug/probe_collision_budget.py`. The probe splits #452's craton
class in two:
- **craton:** a cratonic retreat candidate facing non-cratonic crust;
- **craton_both:** craton on both sides.

It also records each pair's relative rotation |ω_a − ω_b| every step.

    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed <S> --step 2e6 \
        --out <runs>/step2-seed<S>.jsonl
    backend/.venv/bin/python bin/debug/probe_collision_budget.py \
        --summarize <runs>/step2-seed*.jsonl --pool <runs>/step2-seed*.jsonl --json <out>

The four runs reproduce #452's 2 Myr merge counts (29, 24, 18 and 26), so the extended probe
doesn't change the run.

## Results

Share of the double-counted volume's net growth, pooled:

| Class | All pairs | Continental pairs | 30 Myr before merges |
|---|---|---|---|
| Craton on both sides | 73% | 84% | 97% |
| Craton indenting weaker crust | 1% | 1% | −1% |
| Cap | −0.5% | −0.6% | −0.1% |
| Deformable | 27% | 16% | 5% |

Stall force where each class dominates (closing pair-steps):

| Class | Pair-steps | p10 / p50 / p90 (10¹² N/m) | ≤ 10¹³ N/m |
|---|---|---|---|
| Craton on both sides | 6,035 | 0.43 / 2.0 / 7.4 | 94% |
| Craton indenting weaker crust | 105 | 0.56 / 11.6 / 58 | 50% |
| Cap | 93 | 0.40 / 1.4 / 6.8 | 91% |

**Rotation over jams.** The probe followed 493 jammed episodes of at least 3 steps. Relative
rotation started at a median 0.22°/Myr. By the end of the episode it was a median 95% of its
starting value (p10 34%, p90 242%):
- **Halved in only 17%** of episodes.
- **Grew in 45%.**
- **Pairs that later merged** ended at a median 85% of their starting rotation.

## Verdict

- **Craton on both sides makes the blocked growth.** RFC 0011's yield strength scales with
  the lesser of the two sides' craton strengths. Indenters are rare and need more force to
  hold; they keep the same rule.
- **Jammed pairs keep their relative motion.** The contact drag damps closing, but nothing
  brings rotations together. RFC 0011's contact solve gets a friction-bounded tangential
  term, so that merging on matched rotation is reachable.

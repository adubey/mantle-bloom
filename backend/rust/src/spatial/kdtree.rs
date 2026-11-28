use std::cmp::Ordering;

const LEAF_SIZE: usize = 16;
const NO_CHILD: u32 = u32::MAX;
// Median splits keep the depth near log2(n / LEAF_SIZE), so 64 slots cover any u32-sized cloud.
const STACK_DEPTH: usize = 64;

/// How a squared distance is rounded.
///
/// `scipy.spatial.cKDTree` sums `dx², dy², dz²` in that order, but whether each product is
/// rounded before the add depends on how its wheel was compiled: clang's default
/// `-ffp-contract=on` fuses `s += d * d` into a multiply-add on arm64, while baseline x86-64
/// builds round every product. Using the same rounding as the installed scipy keeps every
/// distance, and every threshold and nearest-node choice built on one, bit-identical to the
/// Python fallback (#1042).
pub trait Metric: Copy + Send + Sync + 'static {
    fn d2(a: [f64; 3], b: [f64; 3]) -> f64;
}

/// Every product rounded: `(dx*dx + dy*dy) + dz*dz`.
#[derive(Clone, Copy)]
pub struct Plain;

/// The last two products fused into multiply-adds: `fma(dz, dz, fma(dy, dy, dx*dx))`.
#[derive(Clone, Copy)]
pub struct Fused;

impl Metric for Plain {
    #[inline(always)]
    fn d2(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        dx * dx + dy * dy + dz * dz
    }
}

impl Metric for Fused {
    #[inline(always)]
    fn d2(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        dz.mul_add(dz, dy.mul_add(dy, dx * dx))
    }
}

#[derive(Clone)]
struct KdNode {
    /// Tight bounding box of the node's own points.
    lo: [f64; 3],
    hi: [f64; 3],
    left: u32,
    right: u32,
    start: u32,
    end: u32,
}

/// Balanced static 3D k-d tree over unit-sphere node positions.
///
/// Median splits along each node's widest axis, with every node keeping its points' tight
/// bounding box for pruning: on clustered clouds (land, ocean or a plate's footprint, queried
/// from the whole globe) a box bound prunes cells that a split-plane bound alone cannot,
/// which is why `cKDTree`'s own callers used to need `balanced_tree=False`.
///
/// Exact squared-distance ties break toward the smaller original point index, whatever the
/// tree's layout. `cKDTree` keeps whichever tied point its traversal reaches first, which
/// depends on its leaf size and build options, so the two agree on every query without an
/// exact tie (#1042).
pub struct KdTree {
    nodes: Vec<KdNode>,
    pts: Vec<[f64; 3]>,
    ids: Vec<u32>,
}

#[inline(always)]
fn better(a: (f64, u32), b: (f64, u32)) -> bool {
    a.0 < b.0 || (a.0 == b.0 && a.1 < b.1)
}

fn bounds(items: &[(u32, [f64; 3])]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &(_, p) in items {
        for a in 0..3 {
            if p[a] < lo[a] {
                lo[a] = p[a];
            }
            if p[a] > hi[a] {
                hi[a] = p[a];
            }
        }
    }
    (lo, hi)
}

impl KdTree {
    pub fn build(points: &[[f64; 3]]) -> Self {
        let n = points.len();
        let mut items: Vec<(u32, [f64; 3])> = points
            .iter()
            .enumerate()
            .map(|(i, &p)| (i as u32, p))
            .collect();
        let mut nodes = Vec::with_capacity((n / LEAF_SIZE).max(1) * 2 + 1);
        if n > 0 {
            Self::build_node(&mut items, 0, n, &mut nodes);
        }
        let mut pts = Vec::with_capacity(n);
        let mut ids = Vec::with_capacity(n);
        for (id, p) in items {
            ids.push(id);
            pts.push(p);
        }
        Self { nodes, pts, ids }
    }

    pub fn len(&self) -> usize {
        self.pts.len()
    }

    fn build_node(
        items: &mut [(u32, [f64; 3])],
        start: usize,
        end: usize,
        nodes: &mut Vec<KdNode>,
    ) -> u32 {
        let idx = nodes.len() as u32;
        let (lo, hi) = bounds(&items[start..end]);
        nodes.push(KdNode {
            lo,
            hi,
            left: NO_CHILD,
            right: NO_CHILD,
            start: start as u32,
            end: end as u32,
        });
        if end - start <= LEAF_SIZE {
            return idx;
        }

        let span = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        let axis = if span[1] > span[0] && span[1] >= span[2] {
            1
        } else if span[2] > span[0] && span[2] > span[1] {
            2
        } else {
            0
        };

        let mid = start + (end - start) / 2;
        items[start..end].select_nth_unstable_by(mid - start, |a, b| {
            a.1[axis]
                .partial_cmp(&b.1[axis])
                .unwrap_or(Ordering::Equal)
        });
        let left = Self::build_node(items, start, mid, nodes);
        let right = Self::build_node(items, mid, end, nodes);
        nodes[idx as usize].left = left;
        nodes[idx as usize].right = right;
        idx
    }

    /// Lower bound on `M::d2(p, q)` for every point `q` in node `ni`: `M::d2` of `p`'s per-axis
    /// offsets from the node's box (zero along an axis `p` lies within). For any `q` in the box,
    /// each offset is at most `|p[a] - q[a]|` as rounded, since floating-point subtraction is
    /// monotonic, and squaring and the rounded sums `M` uses are monotonic too, so the bound
    /// never exceeds a point's own computed distance.
    #[inline(always)]
    fn box_bound<M: Metric>(&self, ni: u32, p: [f64; 3]) -> f64 {
        let node = &self.nodes[ni as usize];
        let mut off = [0.0f64; 3];
        for a in 0..3 {
            if p[a] < node.lo[a] {
                off[a] = node.lo[a] - p[a];
            } else if p[a] > node.hi[a] {
                off[a] = p[a] - node.hi[a];
            }
        }
        M::d2(off, [0.0; 3])
    }

    /// Visit every leaf that may hold a point within `limit(state)` squared distance of `p`,
    /// nearer child first, skipping a node whose box bound exceeds the limit both when it is
    /// pushed and again when it is popped (the limit shrinks as neighbours are found). The
    /// bound never exceeds a point's own distance (see `box_bound`), so no point the limit
    /// admits is ever pruned, exact ties included.
    #[inline(always)]
    fn walk<M: Metric, S>(
        &self,
        p: [f64; 3],
        state: &mut S,
        limit: impl Fn(&S) -> f64,
        mut leaf: impl FnMut(&mut S, &[[f64; 3]], &[u32]),
    ) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = [(0u32, 0.0f64); STACK_DEPTH];
        stack[0] = (0, self.box_bound::<M>(0, p));
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let (ni, bound) = stack[sp];
            if bound > limit(state) {
                continue;
            }
            let node = &self.nodes[ni as usize];
            if node.left == NO_CHILD {
                let (s, e) = (node.start as usize, node.end as usize);
                leaf(state, &self.pts[s..e], &self.ids[s..e]);
                continue;
            }
            let (left, right) = (node.left, node.right);
            let (bl, br) = (self.box_bound::<M>(left, p), self.box_bound::<M>(right, p));
            let ((near, bn), (far, bf)) = if bl <= br {
                ((left, bl), (right, br))
            } else {
                ((right, br), (left, bl))
            };
            let lim = limit(state);
            if bf <= lim {
                stack[sp] = (far, bf);
                sp += 1;
            }
            if bn <= lim {
                stack[sp] = (near, bn);
                sp += 1;
            }
        }
    }

    /// The nearest point with squared distance strictly below `ub2`, as `(d2, id)`.
    pub fn nearest<M: Metric>(&self, p: [f64; 3], ub2: f64) -> Option<(f64, u32)> {
        let mut best = (f64::INFINITY, u32::MAX);
        self.walk::<M, _>(
            p,
            &mut best,
            |b| b.0.min(ub2),
            |b, pts, ids| {
                for (q, &id) in pts.iter().zip(ids) {
                    let d2 = M::d2(p, *q);
                    if d2 < ub2 && better((d2, id), *b) {
                        *b = (d2, id);
                    }
                }
            },
        );
        (best.1 != u32::MAX).then_some(best)
    }

    /// Up to `k` nearest points with squared distance strictly below `ub2`, appended to `out`
    /// sorted ascending by `(d2, id)`.
    pub fn knn<M: Metric>(&self, p: [f64; 3], k: usize, ub2: f64, out: &mut Vec<(f64, u32)>) {
        out.clear();
        if k == 0 {
            return;
        }
        self.walk::<M, _>(
            p,
            out,
            |b| if b.len() == k { b[k - 1].0 } else { ub2 },
            |b, pts, ids| {
                for (q, &id) in pts.iter().zip(ids) {
                    let d2 = M::d2(p, *q);
                    if !(d2 < ub2) {
                        continue;
                    }
                    let cand = (d2, id);
                    if b.len() == k {
                        if !better(cand, b[k - 1]) {
                            continue;
                        }
                        b.pop();
                    }
                    let mut pos = b.len();
                    while pos > 0 && better(cand, b[pos - 1]) {
                        pos -= 1;
                    }
                    b.insert(pos, cand);
                }
            },
        );
    }

    /// Ids of every point with squared distance at most `r2` (cKDTree's `query_ball_point`
    /// test), appended to `out` in ascending order.
    pub fn ball<M: Metric>(&self, p: [f64; 3], r2: f64, out: &mut Vec<u32>) {
        let start = out.len();
        self.walk::<M, _>(
            p,
            out,
            |_| r2,
            |b, pts, ids| {
                for (q, &id) in pts.iter().zip(ids) {
                    if M::d2(p, *q) <= r2 {
                        b.push(id);
                    }
                }
            },
        );
        out[start..].sort_unstable();
    }

    /// How many points lie within squared distance `r2`.
    pub fn ball_count<M: Metric>(&self, p: [f64; 3], r2: f64) -> usize {
        let mut count = 0usize;
        self.walk::<M, _>(
            p,
            &mut count,
            |_| r2,
            |c, pts, _| {
                *c += pts.iter().filter(|q| M::d2(p, **q) <= r2).count();
            },
        );
        count
    }

    /// Every point position in original index order.
    pub fn points_by_id(&self) -> Vec<[f64; 3]> {
        let mut out = vec![[0.0; 3]; self.pts.len()];
        for (q, &id) in self.pts.iter().zip(&self.ids) {
            out[id as usize] = *q;
        }
        out
    }

    #[inline]
    pub fn query_nearest(&self, p: [f64; 3]) -> usize {
        self.nearest::<Plain>(p, f64::INFINITY).map_or(0, |(_, id)| id as usize)
    }

    /// Return up to `k` nearest neighbor indices sorted ascending by `(dist2, id)`.
    pub fn query_knn(&self, p: [f64; 3], k: usize) -> Vec<usize> {
        let mut best = Vec::with_capacity(k + 1);
        self.knn::<Plain>(p, k, f64::INFINITY, &mut best);
        best.into_iter().map(|(_, id)| id as usize).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sort_by_distance_then_id(best: &mut [(f64, u32)]) {
        best.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal).then_with(|| a.1.cmp(&b.1)));
    }

    fn cloud(n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let lat = -1.2 + 2.4 * ((i * 73 % n) as f64 / n as f64);
                let lon = -3.0 + 6.0 * ((i * 137 % n) as f64 / n as f64);
                [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()]
            })
            .collect()
    }

    fn brute<M: Metric>(pts: &[[f64; 3]], q: [f64; 3]) -> Vec<(f64, u32)> {
        let mut all: Vec<(f64, u32)> = pts.iter().enumerate().map(|(i, p)| (M::d2(q, *p), i as u32)).collect();
        sort_by_distance_then_id(&mut all);
        all
    }

    #[test]
    fn nearest_and_knn_match_bruteforce() {
        let pts = cloud(200);
        let tree = KdTree::build(&pts);
        for q_idx in [0usize, 17, 99, 199] {
            let q = [pts[q_idx][0] + 0.01, pts[q_idx][1] - 0.02, pts[q_idx][2] + 0.005];
            let all = brute::<Plain>(&pts, q);
            assert_eq!(tree.query_nearest(q), all[0].1 as usize);
            let expected_knn: Vec<usize> = all.iter().take(9).map(|x| x.1 as usize).collect();
            assert_eq!(tree.query_knn(q, 9), expected_knn);
        }
    }

    #[test]
    fn bounded_ball_and_count_match_bruteforce() {
        let pts = cloud(500);
        let tree = KdTree::build(&pts);
        for q_idx in [3usize, 250, 499] {
            let q = [pts[q_idx][0] - 0.004, pts[q_idx][1] + 0.003, pts[q_idx][2]];
            let all = brute::<Fused>(&pts, q);
            let ub2 = all[4].0;
            let mut got = Vec::new();
            tree.knn::<Fused>(q, 8, ub2, &mut got);
            let expected: Vec<(f64, u32)> = all.iter().copied().filter(|x| x.0 < ub2).take(8).collect();
            assert_eq!(got, expected);
            assert_eq!(tree.nearest::<Fused>(q, all[0].0), None);

            let r2 = all[20].0;
            let mut ball = Vec::new();
            tree.ball::<Fused>(q, r2, &mut ball);
            let mut expected_ball: Vec<u32> = all.iter().filter(|x| x.0 <= r2).map(|x| x.1).collect();
            expected_ball.sort_unstable();
            assert_eq!(ball, expected_ball);
            assert_eq!(tree.ball_count::<Fused>(q, r2), expected_ball.len());
        }
    }

    #[test]
    fn exact_ties_pick_the_lowest_index() {
        let base = cloud(64);
        let mut pts = base.clone();
        pts.extend(base.iter().rev().copied());
        let tree = KdTree::build(&pts);
        for (i, q) in base.iter().enumerate() {
            let twin = (2 * base.len() - 1 - i) as u32;
            assert_eq!(tree.nearest::<Plain>(*q, f64::INFINITY), Some((0.0, i as u32)));
            let mut two = Vec::new();
            tree.knn::<Plain>(*q, 2, f64::INFINITY, &mut two);
            assert_eq!(two, vec![(0.0, i as u32), (0.0, twin)]);
        }
    }

    #[test]
    fn empty_tree_finds_nothing() {
        let tree = KdTree::build(&[]);
        let mut out = Vec::new();
        assert_eq!(tree.nearest::<Plain>([1.0, 0.0, 0.0], f64::INFINITY), None);
        tree.knn::<Plain>([1.0, 0.0, 0.0], 3, f64::INFINITY, &mut out);
        assert!(out.is_empty());
        assert_eq!(tree.ball_count::<Plain>([1.0, 0.0, 0.0], 4.0), 0);
    }
}

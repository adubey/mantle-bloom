use std::cmp::Ordering;

const LEAF_SIZE: usize = 16;

#[derive(Clone)]
struct KdNode {
    axis: u8,
    split: f64,
    left: u32,
    right: u32,
    start: u32,
    end: u32,
}

/// Balanced static 3D k-d tree over unit-sphere node positions.
///
/// Exact squared-distance ties break toward the smaller original node index,
/// matching `scipy.spatial.cKDTree` on distinct points.
pub struct KdTree {
    nodes: Vec<KdNode>,
    pts: Vec<[f64; 3]>,
    ids: Vec<u32>,
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

    fn build_node(
        items: &mut [(u32, [f64; 3])],
        start: usize,
        end: usize,
        nodes: &mut Vec<KdNode>,
    ) -> u32 {
        let idx = nodes.len() as u32;
        if end - start <= LEAF_SIZE {
            nodes.push(KdNode {
                axis: 0,
                split: 0.0,
                left: u32::MAX,
                right: u32::MAX,
                start: start as u32,
                end: end as u32,
            });
            return idx;
        }

        let mut min_p = [f64::INFINITY; 3];
        let mut max_p = [f64::NEG_INFINITY; 3];
        for &(_, p) in &items[start..end] {
            for a in 0..3 {
                if p[a] < min_p[a] {
                    min_p[a] = p[a];
                }
                if p[a] > max_p[a] {
                    max_p[a] = p[a];
                }
            }
        }
        let span = [
            max_p[0] - min_p[0],
            max_p[1] - min_p[1],
            max_p[2] - min_p[2],
        ];
        let axis = if span[1] > span[0] && span[1] >= span[2] {
            1
        } else if span[2] > span[0] && span[2] > span[1] {
            2
        } else {
            0
        };

        let mid = start + (end - start) / 2;
        let slice = &mut items[start..end];
        let rel_mid = mid - start;
        slice.select_nth_unstable_by(rel_mid, |a, b| {
            a.1[axis]
                .partial_cmp(&b.1[axis])
                .unwrap_or(Ordering::Equal)
        });
        let split = items[mid].1[axis];

        nodes.push(KdNode {
            axis: axis as u8,
            split,
            left: u32::MAX,
            right: u32::MAX,
            start: start as u32,
            end: end as u32,
        });
        let left = Self::build_node(items, start, mid, nodes);
        let right = Self::build_node(items, mid, end, nodes);
        nodes[idx as usize].left = left;
        nodes[idx as usize].right = right;
        idx
    }

    #[inline]
    fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        dx * dx + dy * dy + dz * dz
    }

    #[inline]
    pub fn query_nearest(&self, p: [f64; 3]) -> usize {
        if self.pts.is_empty() {
            return 0;
        }
        let mut best_d2 = f64::INFINITY;
        let mut best_id = u32::MAX;
        let mut stack = [0u32; 64];
        let mut sp = 1usize;
        stack[0] = 0;

        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp] as usize];
            if node.left == u32::MAX {
                let s = node.start as usize;
                let e = node.end as usize;
                for i in s..e {
                    let d2 = Self::dist2(p, self.pts[i]);
                    let id = self.ids[i];
                    if d2 < best_d2 || (d2 == best_d2 && id < best_id) {
                        best_d2 = d2;
                        best_id = id;
                    }
                }
                continue;
            }
            let diff = p[node.axis as usize] - node.split;
            let (near, far) = if diff <= 0.0 {
                (node.left, node.right)
            } else {
                (node.right, node.left)
            };
            if diff * diff <= best_d2 {
                stack[sp] = far;
                sp += 1;
            }
            stack[sp] = near;
            sp += 1;
        }
        best_id as usize
    }

    /// Return up to `k` nearest neighbor indices sorted ascending by `(dist2, id)`.
    pub fn query_knn(&self, p: [f64; 3], k: usize) -> Vec<usize> {
        if self.pts.is_empty() || k == 0 {
            return Vec::new();
        }
        let mut best: Vec<(f64, u32)> = Vec::with_capacity(k + 1);
        let mut stack = [0u32; 64];
        let mut sp = 1usize;
        stack[0] = 0;

        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp] as usize];
            if node.left == u32::MAX {
                let s = node.start as usize;
                let e = node.end as usize;
                for i in s..e {
                    let d2 = Self::dist2(p, self.pts[i]);
                    let id = self.ids[i];
                    if best.len() < k {
                        best.push((d2, id));
                        if best.len() == k {
                            best.sort_by(|a, b| {
                                a.0.partial_cmp(&b.0)
                                    .unwrap_or(Ordering::Equal)
                                    .then_with(|| a.1.cmp(&b.1))
                            });
                        }
                    } else {
                        let worst = best[k - 1];
                        if d2 < worst.0 || (d2 == worst.0 && id < worst.1) {
                            best[k - 1] = (d2, id);
                            let mut pos = k - 1;
                            while pos > 0 {
                                let cur = best[pos];
                                let prev = best[pos - 1];
                                if cur.0 < prev.0 || (cur.0 == prev.0 && cur.1 < prev.1) {
                                    best.swap(pos, pos - 1);
                                    pos -= 1;
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                continue;
            }
            let worst_d2 = if best.len() < k {
                f64::INFINITY
            } else {
                best[k - 1].0
            };
            let diff = p[node.axis as usize] - node.split;
            let (near, far) = if diff <= 0.0 {
                (node.left, node.right)
            } else {
                (node.right, node.left)
            };
            if diff * diff <= worst_d2 {
                stack[sp] = far;
                sp += 1;
            }
            stack[sp] = near;
            sp += 1;
        }
        if best.len() < k {
            best.sort_by(|a, b| {
                a.0.partial_cmp(&b.0)
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| a.1.cmp(&b.1))
            });
        }
        best.into_iter().map(|(_, id)| id as usize).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearest_and_knn_match_bruteforce() {
        let mut pts = Vec::new();
        for i in 0..200 {
            let lat = -1.2 + 2.4 * ((i * 73 % 200) as f64 / 200.0);
            let lon = -3.0 + 6.0 * ((i * 137 % 200) as f64 / 200.0);
            pts.push([
                lat.cos() * lon.cos(),
                lat.cos() * lon.sin(),
                lat.sin(),
            ]);
        }
        let tree = KdTree::build(&pts);
        for q_idx in [0usize, 17, 99, 199] {
            let q = [pts[q_idx][0] + 0.01, pts[q_idx][1] - 0.02, pts[q_idx][2] + 0.005];
            let got = tree.query_nearest(q);
            let expected = pts
                .iter()
                .enumerate()
                .min_by(|(ia, a), (ib, b)| {
                    let da = KdTree::dist2(q, **a);
                    let db = KdTree::dist2(q, **b);
                    da.partial_cmp(&db).unwrap().then_with(|| ia.cmp(ib))
                })
                .unwrap()
                .0;
            assert_eq!(got, expected);

            let knn = tree.query_knn(q, 9);
            let mut all: Vec<(f64, usize)> = pts
                .iter()
                .enumerate()
                .map(|(i, p)| (KdTree::dist2(q, *p), i))
                .collect();
            all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then_with(|| a.1.cmp(&b.1)));
            let expected_knn: Vec<usize> = all.into_iter().take(9).map(|x| x.1).collect();
            assert_eq!(knn, expected_knn);
        }
    }
}

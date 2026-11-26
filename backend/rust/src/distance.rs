const UNSEEN: u32 = u32::MAX;

/// Indexed 4-ary min-heap over node IDs, keyed by `(dist[node], node)` with in-place decrease-key.
pub struct IndexedHeap {
    heap: Vec<u32>,
    pos: Vec<u32>,
}

impl IndexedHeap {
    pub fn new(n: usize) -> Self {
        Self {
            heap: Vec::with_capacity(n / 4 + 16),
            pos: vec![UNSEEN; n],
        }
    }

    pub fn reset(&mut self) {
        for &node in &self.heap {
            self.pos[node as usize] = UNSEEN;
        }
        self.heap.clear();
    }

    #[inline]
    fn less(a: u32, b: u32, dist: &[f64]) -> bool {
        let da = dist[a as usize];
        let db = dist[b as usize];
        da < db || (da == db && a < b)
    }

    #[inline]
    pub fn push_or_decrease(&mut self, node: u32, dist: &[f64]) {
        let p = self.pos[node as usize];
        if p == UNSEEN {
            let idx = self.heap.len();
            self.heap.push(node);
            self.sift_up(idx, node, dist);
        } else {
            self.sift_up(p as usize, node, dist);
        }
    }

    #[inline]
    pub fn pop(&mut self, dist: &[f64]) -> Option<u32> {
        let top = *self.heap.first()?;
        self.pos[top as usize] = UNSEEN;
        let last = self.heap.pop().unwrap();
        if !self.heap.is_empty() {
            self.sift_down(0, last, dist);
        }
        Some(top)
    }

    #[inline]
    fn sift_up(&mut self, mut idx: usize, node: u32, dist: &[f64]) {
        while idx > 0 {
            let parent = (idx - 1) >> 2;
            let p_node = self.heap[parent];
            if Self::less(node, p_node, dist) {
                self.heap[idx] = p_node;
                self.pos[p_node as usize] = idx as u32;
                idx = parent;
            } else {
                break;
            }
        }
        self.heap[idx] = node;
        self.pos[node as usize] = idx as u32;
    }

    #[inline]
    fn sift_down(&mut self, mut idx: usize, node: u32, dist: &[f64]) {
        let len = self.heap.len();
        loop {
            let first_child = (idx << 2) + 1;
            if first_child >= len {
                break;
            }
            let end_child = (first_child + 4).min(len);
            let mut best_child = first_child;
            let mut best_node = self.heap[first_child];
            for c in (first_child + 1)..end_child {
                let c_node = self.heap[c];
                if Self::less(c_node, best_node, dist) {
                    best_child = c;
                    best_node = c_node;
                }
            }
            if Self::less(best_node, node, dist) {
                self.heap[idx] = best_node;
                self.pos[best_node as usize] = idx as u32;
                idx = best_child;
            } else {
                break;
            }
        }
        self.heap[idx] = node;
        self.pos[node as usize] = idx as u32;
    }
}

fn outside_dijkstra(
    n: usize,
    indptr: &[u32],
    indices: &[u32],
    lengths: &[f64],
    eligible: &[bool],
    seeds: &[bool],
    limit_m: Option<f64>,
) -> (Vec<f64>, Vec<i32>) {
    let mut dist = vec![f64::INFINITY; n];
    let mut sources = vec![-1_i32; n];
    let mut heap = IndexedHeap::new(n);

    for i in 0..n {
        if seeds[i] {
            dist[i] = 0.0;
            sources[i] = i as i32;
            heap.push_or_decrease(i as u32, &dist);
        }
    }

    while let Some(u_u32) = heap.pop(&dist) {
        let u = u_u32 as usize;
        let d = dist[u];
        if let Some(lim) = limit_m {
            if d > lim {
                continue;
            }
        }
        if !eligible[u] {
            continue;
        }
        let src = sources[u];
        let start = indptr[u] as usize;
        let end = indptr[u + 1] as usize;
        for e in start..end {
            let v = indices[e] as usize;
            if !eligible[v] {
                continue;
            }
            let cand = d + lengths[e];
            if let Some(lim) = limit_m {
                if cand > lim {
                    continue;
                }
            }
            if cand < dist[v] || (cand == dist[v] && src < sources[v]) {
                dist[v] = cand;
                sources[v] = src;
                heap.push_or_decrease(v as u32, &dist);
            }
        }
    }

    (dist, sources)
}

pub fn outside_distance(
    n: usize,
    indptr: &[u32],
    indices: &[u32],
    lengths: &[f64],
    eligible: &[bool],
    seeds: &[bool],
    members: &[bool],
    limit_m: Option<f64>,
) -> (Vec<f64>, Vec<i32>) {
    if let Some(lim) = limit_m {
        let (d_out, sources) = outside_dijkstra(n, indptr, indices, lengths, eligible, seeds, Some(lim));
        let mut max_reach = f64::NEG_INFINITY;
        let mut any_member = false;
        for i in 0..n {
            if members[i] {
                any_member = true;
                if d_out[i] > max_reach {
                    max_reach = d_out[i];
                }
            }
        }
        if !any_member || max_reach <= 0.5 * lim {
            return (d_out, sources);
        }
    }
    outside_dijkstra(n, indptr, indices, lengths, eligible, seeds, None)
}

pub fn inside_distance_local(
    n: usize,
    indptr: &[u32],
    indices: &[u32],
    lengths: &[f64],
    positive: &[bool],
) -> Vec<f64> {
    let mut allowed = vec![false; n];
    let mut dist = vec![f64::INFINITY; n];
    let mut heap = IndexedHeap::new(n);
    let mut any_boundary = false;

    for i in 0..n {
        if positive[i] {
            allowed[i] = true;
            continue;
        }
        let start = indptr[i] as usize;
        let end = indptr[i + 1] as usize;
        let mut pos_weight_sum = 0.0;
        for e in start..end {
            if positive[indices[e] as usize] {
                pos_weight_sum += lengths[e];
            }
        }
        if pos_weight_sum > 0.0 {
            allowed[i] = true;
            dist[i] = 0.0;
            heap.push_or_decrease(i as u32, &dist);
            any_boundary = true;
        }
    }

    if !any_boundary {
        return dist;
    }

    while let Some(u_u32) = heap.pop(&dist) {
        let u = u_u32 as usize;
        let d = dist[u];
        let start = indptr[u] as usize;
        let end = indptr[u + 1] as usize;
        for e in start..end {
            let v = indices[e] as usize;
            if !allowed[v] {
                continue;
            }
            let cand = d + lengths[e];
            if cand < dist[v] {
                dist[v] = cand;
                heap.push_or_decrease(v as u32, &dist);
            }
        }
    }

    dist
}

pub fn signed_distance(
    n: usize,
    indptr: &[u32],
    indices: &[u32],
    lengths: &[f64],
    eligible: &[bool],
    positive_seeds: &[bool],
    members: &[bool],
    limit_m: Option<f64>,
) -> (Vec<f64>, Vec<i32>) {
    let ((mut d_out, sources), d_in) = rayon::join(
        || outside_distance(n, indptr, indices, lengths, eligible, positive_seeds, members, limit_m),
        || inside_distance_local(n, indptr, indices, lengths, positive_seeds),
    );
    let mut phi = Vec::with_capacity(n);
    for i in 0..n {
        let di = if positive_seeds[i] { d_in[i] } else { 0.0 };
        if positive_seeds[i] {
            d_out[i] = 0.0;
        }
        phi.push(0.5 * (d_out[i] - di));
    }
    (phi, sources)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_signed_distance_chain() {
        let indptr = vec![0, 1, 3, 5, 6];
        let indices = vec![1, 0, 2, 1, 3, 2];
        let lengths = vec![10.0; 6];
        let eligible = vec![true; 4];
        let seeds = vec![true, true, false, false];
        let members = vec![true, true, false, false];
        let (phi, sources) = signed_distance(4, &indptr, &indices, &lengths, &eligible, &seeds, &members, None);
        assert_eq!(phi, vec![-10.0, -5.0, 5.0, 10.0]);
        assert_eq!(sources, vec![0, 1, 1, 1]);
    }
}

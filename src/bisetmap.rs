use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::hash::Hash;

pub struct BisetMap<L: Hash + Clone + Eq, R: Hash + Clone + Eq> {
    left_to_right: HashMap<L, HashSet<R>>,
    right_to_left: HashMap<R, HashSet<L>>,
}

impl<L: Hash + Clone + Debug + Eq, R: Hash + Clone + Debug + Eq> BisetMap<L, R> {
    pub fn new() -> Self {
        BisetMap {
            left_to_right: HashMap::new(),
            right_to_left: HashMap::new(),
        }
    }

    pub fn insert(&mut self, l: L, r: R) {
        self.left_to_right
            .entry(l.clone())
            .or_insert_with(HashSet::new)
            .insert(r.clone());
        self.right_to_left
            .entry(r)
            .or_insert_with(HashSet::new)
            .insert(l.clone());
    }

    pub fn get_left(&self, l: &L) -> Option<&HashSet<R>> {
        self.left_to_right.get(l)
    }

    pub fn get_right(&self, r: &R) -> Option<&HashSet<L>> {
        self.right_to_left.get(r)
    }

    pub fn remove(&mut self, l: &L, r: &R) {
        match self.left_to_right.get_mut(l) {
            Some(set) => {
                set.remove(r);
                if set.is_empty() {
                    self.left_to_right.remove(l);
                }
            }
            None => (),
        };

        match self.right_to_left.get_mut(r) {
            Some(set) => {
                set.remove(l);
                if set.is_empty() {
                    self.right_to_left.remove(r);
                }
            }
            None => (),
        }
    }

    pub fn exists_left(&self, l: &L) -> bool {
        self.left_to_right.contains_key(l)
    }

    pub fn exists_right(&self, r: &R) -> bool {
        self.right_to_left.contains_key(r)
    }

    pub fn flat_iter(&self) -> impl Iterator<Item = (L, R)> + '_ {
        self.left_to_right
            .iter()
            .flat_map(|(l, r)| r.iter().map(|x| (l.clone(), x.clone())))
    }
}

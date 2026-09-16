"""
Shared LCA implementation.

Owner: Teammate 3 — decide the variant (naive / binary lifting /
Euler tour + sparse table) and implement it here so section 03
(slides/03-python-impl.qmd) can import and demo it directly.
"""


class LCA:
    def __init__(self, tree: dict, root):
        """
        tree: adjacency list, e.g. {1: [2, 3], 2: [4, 5], ...}
        root: root node of the tree
        """
        self.tree = tree
        self.root = root
        # TODO: preprocessing (parent pointers, depths, ancestor
        # tables, Euler tour, etc. depending on chosen variant)

    def query(self, u, v):
        """Return the LCA of nodes u and v."""
        raise NotImplementedError


if __name__ == "__main__":
    # Quick manual test using the same tree as slides/02-hand-trace.qmd
    pass

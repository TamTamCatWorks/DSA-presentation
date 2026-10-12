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
    pass

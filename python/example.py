import robin_py

if __name__ == "__main__":
    print("Examples showing usage of robin_py")

    # creating a Graph in robin
    g = robin_py.AdjListGraph()

    for i in range(10):
        g.AddVertex(i)
        g.AddVertex(i+10)
        g.AddEdge(i, i+10)

    # find the corresponding inlier structures
    max_core_indices = robin_py.FindInlierStructure(
        g, robin_py.InlierGraphStructure.MAX_CORE
    )
    max_clique_indices = robin_py.FindInlierStructure(
        g, robin_py.InlierGraphStructure.MAX_CLIQUE
    )




import performance_benchmark as pb

memory_trace_path = pb.core.memory_trace()
adata = pb.data.read_1M_neurons()

N_cells = [
    500_000,
    1_000_000,
    adata.shape[0],
]

for n_cells in N_cells:
    RData = pb.data.downsample(
            adata,
            n=n_cells,
            random_state=617,
            base_path_for_seurat="../data/1M_neurons/for_seurat",
            name="mouse_neurons_downsample.{}cells".format(n_cells),
        )
    print("\nData Verified: {}".format(RData.verify_paths()))
    seurat_outs = pb.methods.seurat(config=RData, n_iters=4)

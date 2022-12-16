

from ..._core._format_outputs import FormattedResults


def format_scanpy_results(
    run,
    save_filename,
    dataset,
    n_cells,
    run_iter,
    include_runner=False,
    return_saved_outs=False,
):

    formatted = FormattedResults(run)
    formatted.format_scanpy_results(
        dataset=dataset,
        n_cells=n_cells,
        run_iter=run_iter,
        include_runner=include_runner,
    )
    if save_filename:
        formatted.save(filename=save_filename)

    if return_saved_outs:
        return formatted.test_load()[0] # obj only
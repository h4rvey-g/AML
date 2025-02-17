library(sccomp)

test_tcell_composition <- function(sc_tcell) {
    # Create directory
    dir.create("results/108.Tcell/distribution", showWarnings = FALSE, recursive = TRUE)

    # Create count data frame for sccomp with proper dataset column
    # tcell_counts <- sc_tcell@meta.data %>%
    #     select(cell_type_dtl, dataset, group) %>%
    #     mutate(
    #         count = 1
    #     ) %>%
    #     group_by(cell_type_dtl, dataset, group) %>%
    #     summarize(count = n(), .groups = "drop")

    # Run composition test with updated parameters
    sc_tcell <- sc_tcell %>%
        AddMetaData(
            object = .,
            metadata = .$dataset,
            col.name = "sample"
        ) %>%
        AddMetaData(
            object = .,
            metadata = .$cell_type_dtl,
            col.name = "cell_group"
        ) 
    composition_test <- sc_tcell %>%
        sccomp_estimate(
            formula_composition = ~ group,
            # formula_variability = ~1,
            .sample = sample,
            .cell_group = cell_group,
            bimodal_mean_variability_association = TRUE,
            # .abundance = count, # Updated from .count to .abundance
            cores = 20,
            verbose = TRUE
        ) %>%
        sccomp_remove_outliers(cores = 31, verbose = FALSE) %>% 
        sccomp_test()


    # Generate plots
    p1 <- composition_test %>%
        sccomp_boxplot(factor = "group") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/boxplot.png", p1, width = 10, height = 8)

    p2 <- composition_test %>%
        plot_1D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/effect_size.png", p2, width = 8, height = 6)

    p3 <- composition_test %>%
        plot_2D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/abundance_variability.png", p3, width = 8, height = 6)

    # Save fold changes
    fold_changes <- composition_test %>%
        sccomp_proportional_fold_change(
            formula_composition = ~ group,
            from = "normal",
            to = "tumor"
        ) %>%
        select(cell_group, statement)

    write_tsv(fold_changes, "results/108.Tcell/distribution/fold_changes.tsv")

    # Save test results
    write_tsv(
        composition_test %>%
            select(-count_data),
        "results/108.Tcell/distribution/test_results.tsv"
    )

    composition_test
}

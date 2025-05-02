cell_communication <- function(sc_final, sc_mye, sc_tcell) {
    tar_load(sc_final)
    sc_final <- sc_final %>%
        filter(group == "tumor")
    library(liana)
    library(magrittr)
    library(nichenetr)

    # Load myeloid and T cell data if needed
    tar_load(sc_mye)
    tar_load(sc_tcell)

    # Extract cell_type_dtl from myeloid and T cell objects
    mye_meta <- sc_mye@meta.data %>%
        dplyr::select(cell_type_dtl) %>%
        tibble::rownames_to_column("cell_id")

    tcell_meta <- sc_tcell@meta.data %>%
        dplyr::select(cell_type_dtl) %>%
        tibble::rownames_to_column("cell_id")

    # Get current metadata
    current_meta <- sc_final@meta.data %>%
        tibble::rownames_to_column("cell_id")

    # Join the detailed cell types
    updated_meta <- current_meta %>%
        dplyr::left_join(mye_meta, by = "cell_id", suffix = c("", "_mye")) %>%
        dplyr::left_join(tcell_meta, by = "cell_id", suffix = c("", "_tcell"))

    # Update cell_type_dtl where needed
    updated_meta <- updated_meta %>%
        dplyr::mutate(cell_type_dtl = case_when(
            !is.na(cell_type_dtl_tcell) ~ cell_type_dtl_tcell,
            !is.na(cell_type_dtl_mye) ~ cell_type_dtl_mye,
            TRUE ~ cell_type_dtl
        )) %>%
        dplyr::select(-cell_type_dtl_mye, -cell_type_dtl_tcell) %>%
        tibble::column_to_rownames("cell_id")

    # Update sc_final metadata
    sc_final@meta.data <- updated_meta
    liana_test <- liana_wrap(sc_final, idents_col = "cell_type_dtl")
    liana_test <- liana_test %>%
        liana_aggregate()
    library(furrr)
    plan(multisession, workers = 20)

    unique_cell_types <- unique(sc_final$cell_type_dtl)

    liana_dotplot <- function(liana_res,
                              source_groups = NULL,
                              target_groups = NULL,
                              ntop = NULL,
                              specificity = "natmi.edge_specificity",
                              magnitude = "sca.LRscore",
                              y.label = "Interactions (Ligand -> Receptor)",
                              size.label = "Interaction\nSpecificity",
                              colour.label = "Expression\nMagnitude",
                              show_complex = TRUE,
                              size_range = c(2, 10),
                              invert_specificity = FALSE,
                              invert_magnitude = FALSE,
                              facet_by = "source",
                              invert_function = function(x) -log10(x + 1e-10)) {
        if (show_complex) {
            entities <- c("ligand.complex", "receptor.complex")
        } else {
            entities <- c("ligand", "receptor")
        }

        # Modify for the plot
        liana_mod <- liana_res %>%
            # Filter to only the cells of interest
            `if`(
                !is.null(source_groups),
                filter(., source %in% source_groups),
                .
            ) %>%
            `if`(
                !is.null(target_groups),
                filter(., target %in% target_groups),
                .
            )


        if (!is.null(ntop)) {
            # Subset to the X top interactions
            top_int <- liana_mod %>%
                distinct_at(entities) %>%
                head(ntop)
            liana_mod %<>% inner_join(top_int, by = entities)
        }

        if (invert_magnitude) {
            liana_mod %<>% mutate(!!magnitude := invert_function(.data[[magnitude]]))
        }
        if (invert_specificity) {
            liana_mod %<>% mutate(!!specificity := invert_function(.data[[specificity]]))
        }

        liana_mod %<>%
            dplyr::rename(magnitude = !!magnitude) %>%
            dplyr::rename(specificity = !!specificity) %>%
            unite(entities, col = "interaction", sep = " -> ") %>%
            unite(c("source", "target"), col = "source_target", remove = FALSE)



        # ensure levels & order is kept the plot
        interactions_order <- liana_mod %>%
            pull("interaction") %>%
            unique()
        liana_mod %<>%
            mutate(interaction = factor(interaction, levels = rev(interactions_order))) %>%
            mutate(across(where(is.character), as.factor))

        # colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
        cbPalette <- c(
            "#E69F00", "#56B4E9",
            "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#DF69A7"
        )

        # plot
        suppressWarnings(
            ggplot(
                liana_mod,
                aes(
                    x = if (facet_by == "source") target else source,
                    y = interaction,
                    colour = magnitude,
                    size = specificity,
                    group = if (facet_by == "source") target else source
                )
            ) +
                geom_point() +
                scale_color_gradientn(colours = viridis::viridis(20)) +
                scale_size_continuous(range = size_range) +
                facet_grid(. ~ get(facet_by),
                    space = "free",
                    scales = "free",
                    switch = "y"
                ) +
                labs(
                    y = y.label,
                    colour = colour.label,
                    size = size.label,
                    x = if (facet_by == "source") "Target" else "Source",
                    title = stringr::str_to_title(facet_by)
                ) +
                theme_bw(base_size = 20) +
                theme(
                    legend.text = element_text(size = 16),
                    axis.text.x = element_text(
                        colour =
                            cbPalette[1:length(
                                unique(if (facet_by == "source") liana_mod$source else liana_mod$target)
                            )],
                        face = "bold",
                        size = 23
                    ),
                    axis.title.x = element_text(colour = "gray6"),
                    axis.text.y = element_text(
                        size = 18,
                        vjust = 0.5
                    ),
                    legend.title = element_text(size = 18),
                    panel.spacing = unit(0.1, "lines"),
                    strip.background = element_rect(fill = NA),
                    plot.title = element_text(vjust = 0, hjust = 0.5, colour = "gray6"),
                    strip.text = element_text(size = 24, colour = "gray6") # ,
                    # strip.text.y.left = element_text(angle = 0)
                )
        )
    }

    options(future.globals.maxSize = 2000 * 1024^2)
    future_map(unique_cell_types, function(cell_type) {
        # Replace any slashes with underscores in the cell type name for the filename
        cell_type_safe <- gsub("/", "_", cell_type)
        tryCatch(
            {
                p <- liana_test %>%
                    liana_dotplot(
                        source_groups = c(cell_type),
                        ntop = 20
                    ) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(paste0("results/106.fine_analysis/liana_dotplot_send_", cell_type_safe, ".png"), p, width = 12, height = 10)
            },
            error = function(e) {
                message(paste("Error processing cell type", cell_type, ":", e$message))
            }
        )
        # do receiver
        tryCatch(
            {
                p <- liana_test %>%
                    liana_dotplot(
                        target_groups = c(cell_type),
                        ntop = 20,
                        facet_by = "target"
                    ) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(paste0("results/106.fine_analysis/liana_dotplot_receive_", cell_type_safe, ".png"), p, width = 12, height = 10)
            },
            error = function(e) {
                message(paste("Error processing cell type", cell_type, ":", e$message))
            }
        )
    })
    key_types <- c("Adipo", "Fib", "TREM2_LAM", "APOE_LAM")
    p <- liana_test %>%
        liana_dotplot(
            source_groups = c(key_types),
            target_groups = c(key_types),
            ntop = 20
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("results/106.fine_analysis/liana_dotplot_key.png", p, width = 24, height = 10)
    liana_trunc <- liana_test %>%
        # only keep interactions concordant between methods
        filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

    png("results/106.fine_analysis/liana_heatmap.png", width = 12, height = 10, units = "in", res = 300)
    heat_freq(liana_trunc)
    dev.off()
    liana_trunc_filtered <- liana_test %>%
        filter(aggregate_rank <= 0.01, source %in% key_types, target %in% key_types)
    png("results/106.fine_analysis/liana_heatmap_filtered.png", width = 7, height = 6, units = "in", res = 300)
    heat_freq(liana_trunc_filtered)
    dev.off()
    p <- chord_freq(
        liana_trunc
        # source_groups = c("CD8 T", "NK", "B"),
        # target_groups = c("CD8 T", "NK", "B")
    )
    png("results/106.fine_analysis/liana_chord.png", width = 12, height = 10, units = "in", res = 300)
    p
    dev.off()
    p <- chord_freq(
        liana_trunc_filtered
        # source_groups = c("CD8 T", "NK", "B"),
        # target_groups = c("CD8 T", "NK", "B")
    )
    png("results/106.fine_analysis/liana_chord_filtered.png", width = 7, height = 6, units = "in", res = 300)
    p
    dev.off()
    "data/106.fine_analysis"
}

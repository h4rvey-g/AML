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
        rename(magnitude = !!magnitude) %>%
        rename(specificity = !!specificity) %>%
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
                x = target,
                y = interaction,
                colour = magnitude,
                size = specificity,
                group = target
            )
        ) +
            geom_point() +
            scale_color_gradientn(colours = viridis::viridis(20)) +
            scale_size_continuous(range = size_range) +
            facet_grid(. ~ source,
                space = "free",
                scales = "free",
                switch = "y"
            ) +
            # scale_x_discrete(position = "right") +
            labs(
                y = y.label,
                colour = colour.label,
                size = size.label,
                x = "Target",
                title = "Source"
            ) +
            theme_bw(base_size = 20) +
            theme(
                legend.text = element_text(size = 16),
                axis.text.x = element_text(
                    colour =
                        cbPalette[1:length(
                            unique(liana_mod$source)
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

# ------------------------
# author: "Vera Andersson"
# date: "2023-08-27"
# ------------------------

# ------------------------------------------------------------------------------
# Loading packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(igraph)
library(ggraph)
library(gridExtra)

# ------------------------------------------------------------------------------
# Swedish version of data_to_adjacency
# ------------------------------------------------------------------------------

data_to_adjacency <- function(patient_path){
  #' Creates an adjacency matrix with weights as elements (set to 0 for NAs and on
  #' diagonal) given the "patient_path" file path of the raw data, and a data frame
  #' with the node weights. The adjacency matrix and data frame with node weights are
  #' both returned in a list, in that order. Access the adjacency matrix by calling
  #' "data_to_adjacency(patient_path)[[1]]" with your patient file path, and access
  #' the node weights by calling "data_to_adjacency(patient_path)[[2]]" in the same way.
  #' This version has the Swedish symptom names + short descriptions in the adjacency
  #' matrices, which could be used for plotting networks.

  # making adjacency matrix from data
  patient <- read_csv2(patient_path, show_col_types = FALSE)
  colnames(patient)[1] <- "patient_id" # renaming for simplicity
  colnames(patient)[2] <- "relevance" # renaming for simplicity

  ## removing unnecessary columns and adding translated symptoms
  adjacency <- patient[!is.na(patient$relevance),] %>%
    select(-c(starts_with("x"), `Modify?`, Unkn, starts_with("..."))) %>%
    filter_at(1, all_vars(!is.na(.)))

  ## extracting node weights and creating data frame of node weights
  node_weights <- adjacency %>%
    select(patient_id, `Painful?`) %>%
    tail(-1) %>% # not a weight, just a string "Freq" not needed
    mutate(`Painful?` = as.numeric(`Painful?`)) %>%
    rename("node" = patient_id)

  ## removing symptoms not chosen by patient and columns no longer needed
  adjacency <- adjacency %>%
    filter(relevance != 0) %>% # removing symptoms not chosen by patient (row)
    select_if(function(col) !all(col[1] == 0)) %>% # removing symptoms not chosen by patient (col)
    select(-c(relevance, `Painful?`))

  ## selecting the symptoms chosen by patient
  node_weights <- node_weights%>%
    inner_join(adjacency, by = c("node" = "patient_id")) %>%
    select(node, `Painful?`) %>%
    rename("node_weights" = `Painful?`)

  ## removing column and rows no longer needed and converting to numeric
  adjacency <- adjacency[-1,] %>% # not needed
    select(- patient_id) %>%
    mutate_all(as.numeric) %>%
    replace(is.na(.), 0) %>% # there are no edges if weights NA, so we can replace NA with 0
    as.data.frame()

  rownames(adjacency) <- colnames(adjacency) # set row names

  # directed network where rows indicate starting node and columns end node (data is
  # given as the transpose of this)
  adjacency <- t(adjacency)

  list(as.data.frame(adjacency), node_weights) %>%
    return()
}

# ------------------------------------------------------------------------------
# Pruning method 1: Edge betweenness
# ------------------------------------------------------------------------------

pruning_edge_betweenness <- function(adjacency_node_weights, gamma = 0.95) {
  #' Computes the resulting graph when the method is applied to the network
  #' defined by the adjacency matrix and node weights in the list "adjacency_node_weight"
  #' (defined in this way to be used with the function "data_to_adjacency" above).
  #' The method uses the edge betweenness centrality according to Method I in the thesis.
  #' Can tune the value 0 <= "gamma" <= 1 to determine the level of pruning (0 means no
  #' pruning, 1 means that we have the same number of edges as the number of nodes).
  #' Returns a list of the pruned graph object and the original node weights, in that order.

  # extracting adjacency matrix and node weights
  adjacency <- adjacency_node_weights[[1]]
  node_weights_df <- adjacency_node_weights[[2]]

  graph <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                       mode = "directed",
                                       weighted = TRUE)

  g <- graph # to save original "graph" and update "g" in loop below
  n <- ceiling(gamma * (length(E(g)) - length(V(g)))) # number of edges to prune, # edges = # nodes if gamma = 1
  n_edges <- length(E(g))
  edge_weights_keep <- E(g)$weight

  if (n == 0) {
    # if there are no edges to prune, do nothing
    return(list(g, node_weights_df$node_weights))
  }

  # stopping when # nodes = # edges - n (> because we remove an edge within the loop)
  while (length(E(g)) > (n_edges - n)) {
    # transforming edge weights by incorporating node weights
    weights_df <- as_ids(head_of(g, E(g))) %>% # end nodes of each edge
      data_frame() %>%
      left_join(node_weights_df, by = c("." = "node")) %>% # finding end node weights

      # the end nodes are in the order of the nodes in E(g), so E(g)$weight gives the
      # corresponding weight of the edge connected to each end node, whose weight is
      # defined by node_weights

      # normalizing weights and computing new edge weights (NOTE - we divide by sum of E(graph)$weight, not E(g)$weight)
      mutate(new_edge_weights = node_weights / sum(node_weights) * E(g)$weight / sum(E(graph)$weight))

    # creating new weights and replacing the old weights in the graph object g
    new_edge_weights <- weights_df$new_edge_weights

    # edge betweenness with new weights (use 1 - weights because high weight = "short" path and
    # the function "edge_betweenness" interprets the weights as distances, and weights must be non-neg.)
    edge_betweenness_new_weights <- edge_betweenness(g, weights = 1 - new_edge_weights)

    # creating table of all edges in new_g with corresponding indices
    edge_list <- E(g) %>%
      as_ids() %>%
      as_tibble() %>%
      mutate(ind = seq_along(as_ids(E(g))))

    # saving the indices of edges to keep by ordering edge betweenness scores
    best_edges <- order(edge_betweenness_new_weights, decreasing = TRUE) %>%
      head(- 1) # remove worst edge (could have same edge betweenness)

    removed_edge_ind <- order(edge_betweenness_new_weights, decreasing = TRUE) %>%
      tail(1)

    # new pruned subgraph with fewer edges
    subgraph <- subgraph.edges(g,
                               eids = edge_list$ind[-removed_edge_ind],
                               delete.vertices = FALSE) # update g

    # update starting graph and corresponding weights
    g <- subgraph
    edge_weights_keep <- new_edge_weights[edge_list$ind[-removed_edge_ind]]
  }

  E(g)$weight <- edge_weights_keep # updating weights of g

  return(list(g, node_weights_df$node_weights))
}

# --------------------------------------------------------------------------------------
# Pruning method 2: Same number of edges as above in order of edge weights from raw data
# --------------------------------------------------------------------------------------

pruning_og_edge_w <- function(adjacency_node_weights, gamma = 0.95) {
  #' Computes the resulting graph when the method is applied to the network
  #' defined by the adjacency matrix and node weights in the list "adjacency_node_weight"
  #' (defined in this way to be used with the function "data_to_adjacency" above).
  #' The method only keeps the edges with the largest original edge weights, with the
  #' same number of edges kept as for Method 1 above. Can tune the value of
  #' 0 <= "gamma" <= 1 to determine the level of pruning (0 means no pruning,
  #' 1 means that we have the same number of edges as the number of nodes). Returns
  #' a list of the pruned graph object and the original node weights, in that order.

  # extracting adjacency matrix and node weights
  adjacency <- adjacency_node_weights[[1]]
  node_weights_df <- adjacency_node_weights[[2]]

  graph <- graph_from_adjacency_matrix(as.matrix(adjacency),
                                       mode = "directed",
                                       weighted = TRUE)

  g <- graph # to save original "graph" and update "g" in loop below
  n <- ceiling(gamma * (length(E(g)) - length(V(g)))) # number of edges to prune, # edges = # nodes if gamma = 1
  n_edges <- length(E(g))

  if (n == 0) {
    # if there are no edges to prune, do nothing
    return(list(g, node_weights_df$node_weights))
  }

  # creating table of all edges in new_g with corresponding indices
  edge_list <- E(g) %>%
    as_ids() %>%
    as_tibble() %>%
    mutate(ind = seq_along(as_ids(E(g))))

  # saving the indices of edges to remove by ordering by the original edge weights
  removed_edge_ind <- order(E(g)$weight, decreasing = TRUE) %>%
    tail(n)

  pruned_g <- subgraph.edges(g,
                             eids = edge_list$ind[-removed_edge_ind],
                             delete.vertices = FALSE) # update g

  return(list(pruned_g, node_weights_df$node_weights))
}

# ------------------------------------------------------------------------------
# Functions for plotting networks
# ------------------------------------------------------------------------------

wrap_strings <- function(string_vec, width){
  #' Takes a vector of strings and the width of each line as arguments and wraps
  #' the string. Used in plots of networks.
  map_chr(string_vec, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  }
  )
}

plot_weighted_network <- function(g, title, node_weights = 90,
                                  edge_weights = 0.1, text_size = 2.0){
  #' Function for plotting the graph object "g" with the string "title" as the title
  #' of the plot, the nodes sized according to the numerical vector "node_weights",
  #' the edge widths according to the numerical vector "edge_weights", labels of
  #' edges "edge_label", and text size set by the numerical argument "text_size".
  #' Returns the resulting plot.

  size <- node_weights * 0.2
  V(g)$weight <- node_weights * 0.2 # used for length of edges with arrows

  ## apply the function to wrap the node labels
  V(g)$name <- wrap_strings(V(g)$name, 12)

  # producing plot
  plot <- ggraph(g, layout = "stress") +
    ggtitle(title) +
    geom_edge_arc(color = "gray70",
                  strength = - 0.1,
                  arrow = arrow(angle = 15, # width of arrow head
                                length = unit(0.15, "inches"), # length of arrow head
                                ends = "last",
                                type = "closed"),
                  aes(edge_width = edge_weights,
                      start_cap = circle(node1.weight * 1.11, "pt"), # arrows start at node border
                      end_cap = circle(node2.weight + 5, "pt"), # space between end node and arrowhead
                  ),
                  show.legend = c(edge_width = FALSE)) + # edge ending outside of node
    geom_node_point(color = "indianred1",
                    alpha = 0.8, # lower opacity to see edges behind nodes connected to different nodes
                    show.legend = c(size = FALSE, color = FALSE),
                    size = size) +
    geom_node_text(alpha = 1,
                   aes(label = as_ids(V(g))),
                   repel = FALSE,
                   fontface = "bold",
                   colour = "black",
                   size = text_size
    ) +
    scale_edge_width(range = c(0, 1.5)) + # control size
    theme_graph(base_family = "serif", # otherwise error message about font
                plot_margin = margin(90,90,90,90) # used for pdfs, not for other plots
    ) +
    coord_cartesian(clip = "off")
  return(plot)
}

pruned_to_pdf <- function(data_file, save_as, gamma = 0.95) {
  #' Produces a pdf of the resulting networks for each method and the original
  #' network all on the same page with the original network at the top left. The
  #' pdf is saved as "file_name" and the data files of the original data are given
  #' by the list of strings "data_files". The value of "gamma" is set to 0.95.

  # Producing the plots for each method

  adjacencies <- map(data_file, function(x) {
    adjacency_node_weights <- data_to_adjacency(x)
    adjacency <- adjacency_node_weights[[1]]
    node_weights_df <- adjacency_node_weights[[2]]
    return(list(adjacency, node_weights_df))
    })

  ## original plot
  plot_original <- map(adjacencies, function(x) {
    original_g <- graph_from_adjacency_matrix(as.matrix(x[[1]]),
                                              mode = "directed",
                                              weighted = TRUE)
    plot_weighted_network(original_g,
                          "Originalnätverk",
                          edge_weights = E(original_g)$weight * 0.01)
    })

  # edge betweenness network
  plot_edge_b <- map(adjacencies, function(x) {
    graph_pruning_edge_betweenness <- pruning_edge_betweenness(x, gamma)

    plot_weighted_network(graph_pruning_edge_betweenness[[1]],
                             "Edge betweenness")
    })

  # cutoff original edge weights
  plot_og_edge_w <- map(adjacencies, function(x) {
    graph_pruning_og_edge_w <- pruning_og_edge_w(x, gamma)

    plot_weighted_network(graph_pruning_og_edge_w[[1]],
                             "Givna kantvikter")
    })

  # arranging above plots in desired order
  arranged_plots <- c()
  for(i in 1:length(data_file)) {
    name <- paste("plot", i, sep = "")
    arranged_plots <- c(arranged_plots, assign(name, c(plot_edge_b[i],
                                                       plot_og_edge_w[i],
                                                       plot_original[i])))
    }

  # setting desired layout, being 2 + 2 + 2 (3 rows, 2 columns)
  layout <- rbind(c(3,NA),
                  c(2,1))

  # saving the plots
  ggsave(filename = save_as,
         plot = marrangeGrob(arranged_plots,
                             nrow = 2,
                             ncol = 2,
                             layout_matrix = layout,
                             top = quote(paste("Sida", g, "av", npages, "\n", "Patient-ID:",
                                               patientIDs[g]))),
         width = 15, height = 15)
  }

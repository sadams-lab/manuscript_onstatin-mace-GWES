# Based on method descriptions from https://www.biorxiv.org/content/10.1101/353193v1.full

# 1 Ensemble Methods  ----

# * 1.1 Paired Selection Frequency ----

paired_selection_freq <- function(var1, var2, ranger_obj) {
    var_names <- ranger_obj$forest$independent.variable.names
    n_trees <- final_model$forest$num.trees
    trees <- lapply(ranger_obj$forest$split.varIDs, FUN = function(x) {var_names[x]})
    var1_count <- sum(mapply(trees, FUN = function(x, y) {y %in% x}, y = var1 ))
    var2_count <- sum(mapply(trees, FUN = function(x, y) {y %in% x}, y = var2 ))
    var12_actual_on <- sum(mapply(trees, FUN = function(x, y, z) {y %in% x & z %in% x}, z = var1, y = var2 ))
    var12_prob_on <- (var1_count / n_trees) * (var2_count / n_trees)
    print((var12_actual_on/n_trees) - var12_prob_on)
    bi <- binom.test(c(var12_actual_on, n_trees - var12_actual_on), var12_prob_on, alternative = "greater")
    p <- bi$p.value
    return(p)
}


# * 1.2 Split Asymmetry----

# Calculate the standard error per https://www.biorxiv.org/content/10.1101/353193v1.full
calc_se <- function(nABl, nABr, vABl, vABr) {
    return(
        sqrt(
            (((nABl - 1) * vABl) + ((nABr - 1) * vABr)) / (nABl + nABr - 2)
        )
    )
}

# Calculate the t statistic per https://www.biorxiv.org/content/10.1101/353193v1.full
cal_t <- function(nABl, nABr, mABl, mABr, se) {
    return(
        sqrt(
            (nABl * nABr) / (nABl + nABr)) * ((mABl - mABr) / se)
    )
}

# Given a data frame and two variables, create a single node decision tree
# Where vtop is the trunk, and vsub is on both sides of the tree
# Return the difference in numbers assigned to each side of the tree
one_split <- function(vtop, vsub, df) {
    right <- df[df[,vtop] == 0,]
    left <- df[df[,vtop] == 1,]
    right_right <- right[right[,vsub] == 0,]
    right_left <- right[right[,vsub] == 1,]
    left_right <- left[left[,vsub] == 0,]
    left_left <- left[left[,vsub] == 1,]
    return(
        list(
            "right_slope" = (nrow(right_right) / nrow(right_left)),
            "left_slope" = (nrow(left_right) / nrow(left_left))
        ))
}

split_asymmetry <- function(vtop, vsub, sample_frac, df, n_sub) {
    samples <- lapply(1:n_sub, 
                      FUN = function(x) {
                          one_split(vtop, vsub, df = dplyr::sample_frac(df, sample_frac))
                      })
    
    sums_of_slopes_right <- sum(sapply(samples, FUN = function(x) {x[['right_slope']][[1]]}))
    sums_of_slopes_left <- sum(sapply(samples, FUN = function(x) {x[['left_slope']][[1]]}))
    sums_of_squares_slopes_right <- sum(sapply(samples, FUN = function(x) {x[['right_slope']][[1]]^2}))
    sums_of_squares_slopes_left <- sum(sapply(samples, FUN = function(x) {x[['left_slope']][[1]]^2}))
    mean_right <- mean(sapply(samples, FUN = function(x) {x[['right_slope']][[1]]}))
    mean_left <- mean(sapply(samples, FUN = function(x) {x[['left_slope']][[1]]}))
    var_right <- sd(sapply(samples, FUN = function(x) {x[['right_slope']][[1]]})) ^ 2
    var_left <- sd(sapply(samples, FUN = function(x) {x[['left_slope']][[1]]})) ^ 2
    se <- calc_se(nABl = n_sub, nABr = n_sub, vABl = var_left, vABr = var_right)
    tAB <- cal_t(nABl = n_sub, nABr = n_sub, mABl = mean_left, mABr = mean_right, se = se)
    return(list(
        "p" = 2*pt(-abs(tAB), n_sub * 2 - 2), 
        "mean_left" = mean_left, 
        "mean_right" = mean_right, 
        "se" = se,
        "vtop" = vtop,
        "vsub" = vsub,
        "left" = sapply(samples, FUN = function(x) {x[['left_slope']][[1]]}),
        "right" = sapply(samples, FUN = function(x) {x[['right_slope']][[1]]})
    ))
}


# * 1.3 Selection Asymmetry ----

# Re-draw right-side decision tree based on a ranger tree info table
draw_right_tree <- function(v_index, tree, children) {
    first_right <- get_splits(v_index, tree, children)[2]
    first_level <- get_splits(first_right, tree, children)
    new_tree <- tree[v_index:length(tree)]
    second_level_1 <- get_splits(which(new_tree %in% first_level[1])[1], tree, children)
    second_level_2 <- get_splits(which(new_tree %in% first_level[1])[2], tree, children)
    return(c(first_right, first_level, second_level_1, second_level_2))
}

# Re-draw left-side decision tree based on a ranger tree info table
draw_left_tree <- function(v_index, tree, children) {
    first_left <- get_splits(v_index, tree, children)[1]
    first_level <- get_splits(first_left, tree, children)
    new_tree <- tree[v_index:length(tree)]
    second_level_1 <- get_splits(which(new_tree %in% first_level[1])[1], tree, children)
    second_level_2 <- get_splits(which(new_tree %in% first_level[1])[2], tree, children)
    return(c(first_left, first_level, second_level_1, second_level_2))
}

# Given two variants and a ranger info table, 
# return if the sub variant is in the right, left, both, or neither side
what_side <- function(vtop, vsub, trees, children) {
    vtop_index<- which(trees %in% vtop)[[1]]
    right_tree <- draw_right_tree(vtop_index, trees, children)
    left_tree <- draw_left_tree(vtop_index, trees, children)
    in_right <- vsub %in% right_tree
    in_left <- vsub %in% left_tree
    return(list("top" = vtop, "sub" = vsub, "in_right" = in_right, "in_left" = in_left))
}

# Given two variants and a ranger object (with write.forest = TRUE)
# Return indices and tree info tables where both variants are present
tree_co_occurr <- function(var1, var2, ranger_obj) {
    var_names <- ranger_obj$forest$independent.variable.names
    trees <- lapply(ranger_obj$forest$split.varIDs, FUN = function(x) {var_names[x]})
    var12_actual_on <- mapply(FUN = function(x, y, z) {y %in% x & z %in% x}, z = var1, y = var2, x = trees)
    indexes <- which(var12_actual_on %in% 1)
    trees_on <- trees[indexes]
    trees_children <- ranger_obj$forest$child.nodeIDs[indexes]
    return(list("indexes" = indexes, "trees" = trees_on, "children" = trees_children))
}

# Main function for this method
selection_asymmetry <- function(vtop, vsub, df, model) {
    trees_tog <- tree_co_occurr(vtop, vsub, model) # In what trees do the vars occurr together
    tree_list <- mapply(FUN = function(x, y) {what_side(vtop, vsub, x, y)}, x = trees_tog[[2]], y = trees_tog[[3]], SIMPLIFY = FALSE) # Get trees as data frames
    if (length(tree_list) > 0) {
        sub_in_tree <- sum(sapply(tree_list, FUN = function(x) {ifelse(x$in_right | x$in_left, 1, 0)})) # In what trees is the sub var in either the left or right branch
        #s <- sum(df[,vtop] == 1) / nrow(df)
        s_expected <- c(0.5 * sub_in_tree, 0.5 * sub_in_tree)
        s_actual <- c(sum(sapply(tree_list, FUN = function(x) {ifelse(x$in_right, 1, 0)})),
                           sum(sapply(tree_list, FUN = function(x) {ifelse(x$in_left, 1, 0)})))
        p <- chisq.test(cbind(s_expected + 1, s_actual + 1))$p.value
    } else {
        s_expected <- 0
        s_actual <- 0
        p <- 1
    }
    return(list(
        "s_expected" = s_expected, 
        "s_actual" = s_actual, 
        "p" = p))
    
}

# Misc. Supporting Functions ---------------------------------------------------------

# Given a  variant_row and a treeinfo table, return the rows corresponding to the 
# left and right splits
get_splits <- function(v_index, tree, var_row_children) {
    left_node <- var_row_children[[1]][v_index]
    right_node <- var_row_children[[2]][v_index]
    left_var <- tree[left_node + 1]
    right_var <- tree[right_node + 1]
    return(c(left_var, right_var))
}


# Misc. RF Functions ------------------------------------------------------

# Given the mtry and the number of variants, whats the mtry fraction
p_mtry <- function(mtry, num_vars) {
    return(mtry/num_vars)
}

# Given a maximum depth, whats the maximum
# number of nodes in a given tree
num_nodes <- function(max_depth) {
    return(2 ^ (max_depth) - 1)
}

# Probability that a random, non-associated variable will be anywhere in a tree assuming no association
prob_in_tree <- function(max_depth, mtry, num_vars) {
    nodes <- num_nodes(max_depth) # How many nodes (maximum possible)
    if (nodes > num_nodes(mtry)) {
        return(NaN)
    }
    p <- (nodes/mtry) * p_mtry(mtry, num_vars)
    if (p >= 1.0) {
        return(1.0)
    } else {
        return(p)
    }
}

fisher_p_method <- function(p_vals) {
    x2 <- -2 * sum(log(p_vals))
    return(pchisq(x2, df=2 * length(p_vals), lower.tail = FALSE))
}

# Main Method -------------------------------------------------------------

# This runs all three ensemble methods for a pair of variants

validate_epistasis <- function(v1, v2, model) {
    # #print("Running Regression")
    # linreg <- summary(glm(data = df, as.formula(paste("MACE_EVENT ~ ", paste(v1, v2, sep = "*"), sep = "")) , family = "binomial"))
    # lin_v1v2_b <- linreg$coefficients[nrow(linreg$coefficients),1]
    # lin_v1v2_p <- linreg$coefficients[nrow(linreg$coefficients),4]

    
    #print("Getting Paired Selection Frequency")
    m1_p <- paired_selection_freq(v1, v2, model)
    
    #print("Checking Split Asymmetry")
    #m2_1 <- split_asymmetry(vtop = v1, vsub = v2, sample_frac = 0.5, df = df, n_sub = 10)
    #m2_2 <- split_asymmetry(vtop = v2, vsub = v1, sample_frac = 0.5, df = df, n_sub = 10)
    
    #print("Checking Selection Asymmetry")
    m3_1 <- selection_asymmetry(vtop = v1, vsub = v2, df = df, model = model)
    m3_2 <- selection_asymmetry(vtop = v2, vsub = v1, df = df, model = model)
    
    return(list(
        "m1_1_p" = m1_p,
        #"m2_1" = m2_1,
        #"m2_2" = m2_2,
        "m3_1" = m3_1,
        "m3_2" = m3_2
        #"lin_v1v2_b" = lin_v1v2_b,
        #"lin_v1v2_p" = lin_v1v2_p
    ))
}


plot_upset = function(lt, title = NULL, value_fun = length, remove_full = FALSE) {

    n = length(lt)
    nm = names(lt)
    
    set_mat = matrix(FALSE, nrow = n, ncol = sum(choose(n, 1:n)))
    rownames(set_mat) = nm
    j = 1
    for(k in 1:n) {
        comb = combn(n, k)
        for(i in 1:ncol(comb)) {
            set_mat[comb[, i], j] = TRUE
            j = j + 1
        }
    }
    if(remove_full) {
    	set_mat = set_mat[, -ncol(set_mat)]
    }

    do_set = function(lt, do_intersect = rep(TRUE, length(lt)), value_fun = length) {
        set1_index = which(do_intersect)
        set2_index = which(!do_intersect)

        s = lt[[ set1_index[1] ]]
        
        for(i in set1_index[-1]) {
            s = intersect(s, lt[[ i ]])
        }

        for(i in set2_index) {
            s = setdiff(s, lt[[ i ]])
        }
        value_fun(s)
    }

    set_value = numeric(ncol(set_mat))
    for(i in seq_len(ncol(set_mat))) {
        # qqcat("intersection: @{paste(set_mat[, i]+0, collapse = ', ')}\n")
        set_value[i] = do_set(lt, set_mat[, i], value_fun = value_fun)
    }
    set_name = sapply(seq_len(ncol(set_mat)), function(i) {
        nm2 = nm[set_mat[, i]]
        paste(nm2, collapse = "&")
    })

    make_upset(set_mat, set_value, title = title)

    return(invisible(structure(set_value, names = set_name)))
}


make_upset = function(set_mat, set_value, title = NULL) {
    set_size = sapply(seq_len(nrow(set_mat)), function(i) {
        sum(set_value[set_mat[i, ]])
    })

    nc = ncol(set_mat)
    nr = nrow(set_mat)
    sn = rownames(set_mat)

    grid.newpage()
    sn_width = max_text_width(sn) + unit(1, "cm")
    pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2, width = unit(c(2, 1), "null"), height = unit(c(4, 1), "null")), x = sn_width, just = "left", width = unit(0.95, "npc") - sn_width, height = unit(0.95, "npc")))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, xscale = c(0, nc), yscale = c(0, max(set_value))))
    grid.rect(1:nc - 0.5, 0, width = 0.6, height = set_value, default.units = "native", just = "bottom", gp = gpar(fill = "black"))
    grid.yaxis(gp = gpar(fontsize = 8))
    grid.text("#intersections", unit(-1.5, "cm"), just = "bottom", rot = 90)
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1, xscale = c(0, nc), yscale = c(0, nr)))
    for(i in seq_len(ncol(set_mat))) {
        grid.points(rep(i - 0.5, nr), 1:nr - 0.5, default.units = "native", size = unit(3, "mm"), pch = 16, gp = gpar(col = ifelse(set_mat[, i], "black", "#CCCCCC")))
        if(sum(set_mat[, i]) >= 2) {
            i_min = min(which(set_mat[, i]))
            i_max = max(which(set_mat[, i]))
            grid.lines(c(i - 0.5, i - 0.5), c(i_min, i_max) - 0.5, default.units = "native", gp = gpar(col = "black", lwd = 2))
        }
    }
    grid.text(sn, x = unit(0, "npc") - unit(2, "mm"), y = 1:nr - 0.5, just = "right", default.units = "native")
    popViewport()
    
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    pushViewport(viewport(x = unit(5, "mm"), width = unit(1, "npc") - unit(5, "mm"), just = "left", xscale = c(0, max(set_size)), yscale = c(0, nr)))
    grid.rect(unit(0, "npc"), 1:nr - 0.5, width = set_size, height = 0.6, just = "left", default.units = "native", gp = gpar(fill = "black"))
    grid.xaxis(main = FALSE, gp = gpar(fontsize = 8))
    grid.text("set size", y = unit(1, "npc") + unit(1, "cm"), just = "bottom")
    popViewport()
    popViewport()

    if(!is.null(title)) {
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
        grid.text(title, gp = gpar(fontsize = 16))
        popViewport()
    }
 
    popViewport()
}

## 7/6/2022

## Visium helper function


## functions defined to assist with plotting and analysis of visium data


#### connectedComponentsFunctions ####


# first define knn function then define the connected components function


##### knn_from_coord #####
# define function to create knn graph from coordinate matrix
# input is the cell coordinates
## this is using code from Giotto: https://github.com/RubD/Giotto/blob/master/R/NN_network.R
knn_from_coord = function(matrix_to_use, cell_names, 
                          k, distance_thresh = NULL){
  
  # define cell names, essentially just an index
  # cell_names = 1:nrow(matrix_to_use) 
  # names(cell_names) = 1:nrow(matrix_to_use)
  # this is specified on the original segmented image data
  names(cell_names) = cell_names
  nn_network = dbscan::kNN(x = matrix_to_use, k = k, sort = TRUE, search = 'kdtree')
  
  # data.table variables
  from = to = weight = distance = from_cell_ID = to_cell_ID = shared = NULL
  
  nn_network_dt = data.table::data.table(from = rep(1:nrow(nn_network$id), k),
                                         to = as.vector(nn_network$id),
                                         weight = 1/(1 + as.vector(nn_network$dist)),
                                         distance = as.vector(nn_network$dist))
  # print the median distance between neighbors
 # print(paste("median distance between neighbors is:", nn_network_dt$distance %>% median()))
  ## this is not correct, I don't understand what kind of distance I'm getting here
  
  nn_network_dt[, from_cell_ID := cell_names[from]]
  nn_network_dt[, to_cell_ID := cell_names[to]]
  if(!(is.null(distance_thresh))){
    nn_network_dt = nn_network_dt[distance < distance_thresh]
  }
  ## convert to igraph object
  all_index = unique(x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID))
  
  # need an undirected graph here
  nn_network_igraph = igraph::graph_from_data_frame(nn_network_dt[,.(from_cell_ID, to_cell_ID, weight, distance)], directed = TRUE, vertices = all_index)
  
  return(nn_network_igraph)
}


# testing
# # generate igraph object from segmented cell coordinates using the 6 nearest neighbors
# test.graph = knn_from_coord(GetTissueCoordinates(d$Ctrl1_A1), 
#                             colnames(d$Ctrl1_A1), 6, 7)
# test.comp = igraph::components(test.graph)
# test.comp.df = reshape2::melt(as.matrix(test.comp$membership)) %>% 
#                        dplyr::select(barcode = Var1, component = value)
# 
# coord_dist = as.matrix(dist(GetTissueCoordinates(d$Ctrl1_A1)))
# coord_dist[test.graph %>% as.matrix()] %>% hist()





##### separate_connected_components ####
# returns a vector with spot component memberships with barcodes as names
# define a function that gets connected components of certain class labels
# this will be for defining certain regions in the visium sections, primarily tumor regions and maybe IFN and inflammatory cytokine loci
separate_connected_components = function(sobj, #label, 
                                         n_neighs = 6){
  
  # calculate spatial neighbors, make adjacency matrix
  coords = GetTissueCoordinates(sobj)
  
  # find neighbors
  ## this will be done by taking the 6 nearest neighbors for each 6 neighbor max 
  ## spot and then taking the median of that * 1.3 (as is done in sqiudpy) to 
  ## restrict to only neighbors within direct contact of the spot
  ## this is radius multiplied by 5 because the visium (V1) has a spot 50um in diameter
  ## with 100um spot center to center distance. Therefore a distance of 125 should be a good
  ## threshold for pruning out things that are too far but I can check
  #dist_thresh = SeuratObject::Radius(sobj@images$slice1) * 5
  
  
  
  # # for each group in label category
  # cc_list = lapply(sobj@meta.data[,label] %>% unique(), 
  #                  function(group){
  #                    group_idx = sobj@meta.data[,label] == group
  #                    coords_f = coords[group_idx,]
  #                    knn_graph = knn_from_coord(coords_f, colnames(sobj[,group_idx]), n_neighs)#, distance_thresh = dist_thresh)
  #                    conn_comp = igraph::components(knn_graph)
  #                    return(conn_comp)
  #                  })
  
  #coords_f = coords[group_idx,]
  
  #knn_graph = knn_from_coord(coords, colnames(sobj), n_neighs)#, distance_thresh = dist_thresh)
  # get the mean path length and multiply by 1.2 to get the distance threshold
  #dist_thresh = igraph::mean_distance(knn_graph) * 1.2
  ## this is around what it should be for visium I think, seems to work alright
  dist_thresh = 8
  #print(dist_thresh)
  knn_graph = knn_from_coord(coords, colnames(sobj), n_neighs, distance_thresh = dist_thresh)                
  conn_comp = igraph::components(knn_graph)
  
  # merge to original order of cells and return the component labels
  comp_df = reshape2::melt(as.matrix(conn_comp$membership)) %>% 
    dplyr::select(barcode = Var1, component = value) %>%
    mutate(barcode = as.character(barcode))
  
  comp_ret = dplyr::left_join(data.frame(barcode = colnames(sobj)), 
                              comp_df, by = 'barcode' ) 
  comp_vec = comp_ret$component
  names(comp_vec) = comp_ret$barcode
  return(comp_vec)
  
  #return(conn_comp)
  #return(cc_list)
  #dist_thresh 
}




# function to get borders around each component and plot them so that I can easily highlight the regions

##### groupBorder ####
# function that produces a border around a set of points using the ahull method
#  need to iterate over alpha values to draw borders until a circle is formed
# by connecting in igraph
# input: data frame resulting from the dbscan cluster annotation
# value: will return a list of sp polygon objects for each procode with a defined cluster
## this will also include sizes for the polygons in a named list indexed by polygon or size

groupBorder <- function(dat, group_label = 'component', alpha_start = 80){
  
  # list of groups to draw borders for
  group.list <- unique(dat@meta.data[,group_label])
  #print(group.list)
  
  # loop over each group to draw borders and get out sp polygon objects
  group.pol <- lapply(group.list, function(g){
    #print(paste("group", g))
    # filter for the procode
    dat.f <- dat[,dat@meta.data[,group_label] == g]
    #print(dat.f)
    
    # loop over function creating border and polygon around a set of points
    # until I have an alpha value that works
    # if(nrow(dat.f) > 3){
    #   repeat{
    #     print(paste("drawing borders w/ alpha", alpha_start))
    #     pol <- try(borderPolygon(dat.f, alpha_start))
    #     
    #     if(!(inherits(pol, "try-error")))
    #       break
    #     
    #     if(alpha_start > 300)
    #       break
    #     
    #     alpha_start <- alpha_start + 10
    #   }
    # } else{ 
    #   
    #   pol = st_convex_hull(st_multipoint(Seurat::GetTissueCoordinates(dat.f) %>% as.matrix()))
    #   
    # }
    
    #st_obj = st_convex_hull(st_multipoint(Seurat::GetTissueCoordinates(dat.f) %>% as.matrix()))
    
    coord = Seurat::GetTissueCoordinates(dat.f, image = 'slice1') %>% as.matrix()
    pol = st_convex_hull(st_multipoint(coord[,c(2,1)]))
    # there is a transform that I have to do for Seurat plotting that I am adding in here
    #coord[,1] = max(GetTissueCoordinates(d$DT1_A1)[,'imagerow']) - imagecol + min(GetTissueCoordinates(d$DT1_A1)[,'imagerow'])
    #pol = st_multipoint(coord[,c(2,1)])
    #pol = st_multipoint(Seurat::GetTissueCoordinates(dat.f, image = 'slice1') %>% as.matrix())
    
    #pol = st_boundary(st_multipoint(Seurat::GetTissueCoordinates(dat.f) %>% as.matrix()))
    #pol = st_polygon(Seurat::GetTissueCoordinates(dat.f) %>% as.matrix())
    # convert to lined polygon or something
    #st_obj = st_segmentize(st_obj, Inf)
    #pol = st_buffer(st_obj, 1.1)
    ## replacing above code with st_buffer to see if that is much easier
    # get the size of the created polygon
    r <- list(polygon = pol, 
              size = st_area(pol))
    return(r)
  }
  )
  names(group.pol) <- group.list
  return(group.pol)
  
}



##### borderPolygon ####

# function to draw alpha shape around points and return polygon vector
# for now I will just get the lines back so I can plot and see if things
# are working

borderPolygon <- function(dat, alpha_start){
  
  # draw the alpha shape around this set of points
  coords = Seurat::GetTissueCoordinates(dat) %>% as.matrix()
  pc.shape <- alphahull::ashape(coords, 
                                alpha = alpha_start)
  # assemble the points into connected lines
  pc.g <- graph.edgelist(cbind(as.character(pc.shape$edges[, "ind1"]), 
                               as.character(pc.shape$edges[, "ind2"])), directed = FALSE)
  
  # run tests to make sure alpha value used forms the correct shape of graph
  ## should change this to return things based on what happened to direct the 
  ## above function to iterate again until it gets the correct value
  ### for now everything will just stop if one of them is messed up
  if (!is.connected(pc.g)) {
    stop("Graph not connected")
  }
  if (any(degree(pc.g) != 2)) {
    stop("Graph not circular")
  }
  if (clusters(pc.g)$no > 1) {
    stop("Graph composed of more than one circle")
  }
  
  # I will now just return the xy coords of the polygon for plotting in ggplot
  cutg = pc.g - E(pc.g)[1]
  # find chain end points
  ends = names(which(degree(cutg) == 1))
  path = unlist(get.shortest.paths(cutg, ends[1], ends[2])[[1]])
  # this is an index into the points
  pathX = as.numeric(V(pc.g)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  
  # # get the xy coords of the points for drawing a polygon
  #pc.poly <- pc.shape$x[pathX, ]
  # # change the colnames to be compatible with my plotting
  #colnames(pc.poly) <- c('Centroid X px', 'Centroid Y px')
  
  # now for this I will return an sf polygon object because this will interface
  # with other spatial testing procedures
  # for plotting this can always be converted back to a matrix by using 
  # the st_coordinates function
  pc.poly <- st_polygon(list(pc.shape$x[pathX,]))
  return(pc.poly)
  
}


##### SpatialFeaturePolyContext ####
# this function plots the defined tumor lesions in the context of the visium image
# can specify a vector of colors for the polygons as border color and this will color each one individually based on that
SpatialFeaturePolyContext = function(d, feature, polys, min_cut, max_cut, border_color = 'white', line_size = .5){
  
  gp = SpatialFeaturePlot(d, feature, image.alpha = 0,min.cutoff = min_cut, max.cutoff = max_cut) + scale_fill_viridis_c(option = "B")
  
  if(length(border_color) == 1){
    border_color_idx = rep(border_color, length(polys))
  } else{
    border_color_idx = border_color
  }
  # add each tumor boundary
  for(i in seq(from = 1, to = length(polys), by = 1)){
    #print(polys[[i]]$polygon)
    gp = gp + geom_polygon(st_coordinates(polys[[i]]$polygon) %>%
                             as.data.frame() %>% dplyr::rename(imagerow = X,
                                                               imagecol = Y) %>% 
                             mutate(!! feature := rep(1, n()), 
                                    # the imagecol vs imagerow gets flipped around here but this ends up working with what I have
                                    # I flip the order of the columns in the coordinate matrix so this ends up being what is needed
                                    imagecol = max(GetTissueCoordinates(d)[,'imagerow']) - imagecol + min(GetTissueCoordinates(d)[,'imagerow'])
                                    #imagecol = max(GetTissueCoordinates(d$DT1_A1)[,'imagecol']) - imagecol + min(GetTissueCoordinates(d$DT1_A1)[,'imagecol'])
                             ),
                           mapping=aes(x = imagerow, y = imagecol),
                           linetype = 2, fill = NA, color = border_color_idx[as.numeric(i)], size = line_size )
  }
  gp
}


##### SpatialFeaturePlotListPoly ####
# define a function that will get the plotting right for the factors with the viridis colors I want
SpatialFeaturePlotListPoly = function(section, features, feature_max, polys, poly_color = NULL){
  
  gplots = lapply(features, function(x){ 
    if(is.null(poly_color)){
      SpatialFeaturePolyContext(section, x, polys, min_cut = -0.1, 
                                max_cut = feature_max[x]
      ) + scale_fill_viridis_c(option = "B")  
    } else {
      SpatialFeaturePolyContext(section, x, polys, min_cut = -0.1, 
                                max_cut = feature_max[x], border_color = poly_color
      ) + scale_fill_viridis_c(option = "B")  
    }
    
    
  })
  
  n <- length(gplots)
  nCol <- floor(sqrt(n))
  grobs <- lapply(gplots, ggplotGrob)
  do.call("grid.arrange", c(grobs, ncol = 3))
  
  
  
}



#### basic plotting ####

##### SpatialFeaturePlotList ####
# define a function that will get the plotting right for the factors with the viridis colors I want
SpatialFeaturePlotList = function(section, features, feature_max){
  
  gplots = lapply(features, function(x){ 
    SpatialFeaturePlot(section, x, image.alpha = 0, min.cutoff = -0.1, 
                       max.cutoff = feature_max[x] 
    ) + scale_fill_viridis_c(option = "B")
    
  })
  
  n <- length(gplots)
  nCol <- floor(sqrt(n))
  grobs <- lapply(gplots, ggplotGrob)
  do.call("grid.arrange", c(grobs, ncol = 3))
  
  
  
}


##### SpatialFeaturePlotListTissue ####
# define a function that will get the plotting right for the factors with the viridis colors I want
# show tissue as well
# SpatialFeaturePlotListTissue = function(section, features, feature_max){
#   
#   gplots = lapply(features, function(x){ 
#     SpatialFeaturePlot(subset(section, cells = colnames(section)[section@assays$SSM], x, image.alpha = 1, min.cutoff = -0.1, 
#                        max.cutoff = feature_max[x], pt.size.factor = 1.1, alpha = .4
#     ) + scale_fill_viridis_c(option = "B")
#     
#   })
#   )
#   n <- length(gplots)
#   nCol <- floor(sqrt(n))
#   grobs <- lapply(gplots, ggplotGrob)
#   do.call("grid.arrange", c(grobs, ncol = 3))
#   
#   
#   
# }

##### SpatialFeaturePlotListTissue ####
# define a function that will get the plotting right for the factors with the viridis colors I want
# show tissue as well
# only positive spots shown
SpatialFeaturePlotListTissuePos = function(section, features, feature_max){
  
  gplots = lapply(features, function(x){ 
    # check how many cells are positive. Must be greater than 3 for plotting
    if(sum(section@assays$SSM@counts[x,]) >= 3){
      SpatialFeaturePlot(subset(section, cells = colnames(section)[section@assays$SSM@counts[x,] %>% as.logical()]), x, image.alpha = 1, min.cutoff = -0.1, 
                         max.cutoff = feature_max[x], pt.size.factor = 1, alpha = .5
      ) + scale_fill_viridis_c(option = "B", limits = c(0,1))      
    } else{
      SpatialFeaturePlot(section, x, image.alpha = 1, min.cutoff = -0.1, 
                         max.cutoff = feature_max[x], pt.size.factor = 1.1, alpha = .5
      ) + scale_fill_viridis_c(option = "B", limits = c(0,1)) 
    }

    
  })
    #gplots = gplots[!(is.null(gplots))]
    n <- length(gplots)
    nCol <- floor(sqrt(n))
    grobs <- lapply(gplots, ggplotGrob)
    do.call("grid.arrange", c(grobs, ncol = 3))
    
    
    
}

# define a function that will get the plotting right for the factors with the viridis colors I want
SpatialPlotSections = function(sections, feature, color_palette){
  
  gplots = lapply(sections, function(x){ 
    SpatialPlot(x, group.by = feature, image.alpha = 0#, #min.cutoff = -0.1, 
                #max.cutoff = feature_max[x] 
    ) + color_palette
    
  })
  
  n <- length(gplots)
  nCol <- floor(sqrt(n))
  grobs <- lapply(gplots, ggplotGrob)
  do.call("grid.arrange", c(grobs, ncol = 4))
  
  
  
}

# define a function that will get the plotting right for the factors with the viridis colors I want
SpatialFeaturePlotSections = function(sections, feature, feature_max, feature_min, color_palette){
  
  gplots = lapply(sections, function(x){ 
    SpatialFeaturePlot(x,  feature, image.alpha = 0, min.cutoff = feature_min, 
                       max.cutoff = feature_max 
    ) + color_palette#scale_fill_viridis_c(option = "B")
    
  })
  
  n <- length(gplots)
  nCol <- floor(sqrt(n))
  grobs <- lapply(gplots, ggplotGrob)
  do.call("grid.arrange", c(grobs, ncol = 4))
  
  
  
}

#### SpatialFeaturePlotSectionsPoly ####
# function to plot several sections with boundaries at the same time for a given feature

SpatialFeaturePlotSectionsPoly = function(sections, feature, feature_max, feature_min, color_palette, poly){
  
  #feature, polys, min_cut, max_cut, border_color = 'white', line_size = .5
  gplots = lapply(names(sections), function(x){ 
    #print(str(poly[x]))
    SpatialFeaturePolyContext(sections[[x]],  feature, poly[[x]], min_cut = feature_min, 
                       max_cut = feature_max 
    )  
    
  })
  
  n <- length(gplots)
  nCol <- floor(sqrt(n))
  grobs <- lapply(gplots, ggplotGrob)
  do.call("grid.arrange", c(grobs, ncol = 4))
  
  
  
}


# Function that will plot tumor lesions identified in metadata by the label_key and label them in their centroid
## I need to finish this, the coordinates are not getting translated properly and it is messing up the plotting
SpatialDimPlotTumLabel = function(section, label_key){
  
  tumor_lesions = unique(section@meta.data[,label_key], )
  tumor_lesions = tumor_lesions[!(is.na(tumor_lesions))]
  # get tumor centroids
  coord = lapply(tumor_lesions, function(x){
    Seurat::GetTissueCoordinates(section[,section@meta.data[,label_key] == x], image = 'slice1') %>% as.matrix()
  })
  names(coord) = tumor_lesions
  
  #print("coords are")
  #print(coord)
  cent = do.call(rbind, lapply(tumor_lesions, function(x){
    #print(paste("tumor_lesion is ", x))
    #print(paste("coords are:", coord[[x]]))
    cret = colMeans(coord[[x]]) %>% t() %>% as.data.frame()
    colnames(cret) = c('imagerow', 'imagecol')
    cret$tumor_lesion = x
    cret
  }
  ))
  print(head(cent))
  
  max_imagecol = max(GetTissueCoordinates(section)[,'imagecol'])
  min_imagecol = min(GetTissueCoordinates(section)[,'imagecol'])
  
  # create df of transformed coordinates with labels
  ## need to set ident to 1 bc of the plotting
  cent.t = mutate(cent, imagerow = max_imagecol - imagerow + min_imagecol, 
                  ident = rep('1', n()))
  
  g = SpatialDimPlot(section)
  
  # add labels
  for(i in tumor_lesions){
    g = g + geom_text(data = cent.t, 
                      mapping=aes(x = imagerow, y = imagecol, label = tumor_lesion),
                      color = 'white')
  }
  g
  
  
}

# get xy coordinates that will match what seurat has for plotting on top of their plot
transform_coordinates_seurat_image = function(section){
  ## need to finish this ##
  imagecol = max(GetTissueCoordinates(section)[,'imagerow']) - imagecol + min(GetTissueCoordinates(section)[,'imagerow'])
}

#### spatialFeaturePlotMulticolor ####
# function that will plot multiple colors at different intensities
SpatialFeaturePlotList = function(section, features, feature_max){
  
  gplots = lapply(features, function(x){ 
    SpatialFeaturePlot(section, x, image.alpha = 0, min.cutoff = -0.1, 
                       max.cutoff = feature_max[x] 
    ) + scale_fill_viridis_c(option = "B")
    
  })
  
  n <- length(gplots)
  nCol <- floor(sqrt(n))
  grobs <- lapply(gplots, ggplotGrob)
  do.call("grid.arrange", c(grobs, ncol = 3))
  
  
  
}


#### runBasicLimmaDE ####
runBasicLimmaDE <- function(exp, mod, cutoff = F){
  
  if(cutoff != F){
    print(paste("filtering expression values on cutoff of", cutoff))
    #exp <- filterExpressedAQ(exp, cutoff = cutoff)
  }
  
  design <- model.matrix(~ mod)
  fit <- lmFit(exp, design)
  fit <- eBayes(fit)
  res <- topTable(fit, n = nrow(exp))
  sample_averages <- t(apply(exp[rownames(res),], 1, function(x) {
    tapply(x, mod, mean)
  }))
  res <- cbind(res, sample_averages)
  res$gene <- rownames(res)
  return(res)
}






#### BayesPrism helpers ####

##### sum_genes #####

# function to sum the expression of a set of genes in each spot
# input a list and returned will be the sum of genes in each list entry
# named by the list entry name
# input gene list and the visium data
# need to make sure the visium data is unnormalized counts

sum_gene_in_spots = function(section, gene_list, slot = 'counts'){
  if(slot == 'counts') dat = section@assays$Spatial@counts %>% as.matrix()
  if(slot == 'data') dat = section@assays$Spatial@data %>% as.matrix()
  
  sum_array = matrix(nrow = ncol(dat), 
                     ncol = length(gene_list), 
                     dimnames = list(colnames(dat),
                                     names(gene_list)))
  
  lapply(names(gene_list), function(idx){
    
    sum_array[,idx] <<- apply(dat[rownames(dat) %in% gene_list[[idx]],], 2, sum)
    
  })
  # return the sum array
  return(sum_array)
  
}  

##### mean_zscores_spots #####

# function to get the mean of Z-scores for a gene list
# input a list and returned will be the sum of genes in each list entry
# named by the list entry name
# input gene list and the visium data
## will update this function if I need it but I can also try module scores

# mean_zscores_spots = function(section, gene_list, slot = 'counts'){
#   if(slot == 'counts') dat = section@assays$Spatial@counts %>% as.matrix()
#   if(slot == 'data') dat = section@assays$Spatial@data %>% as.matrix()
#   
#   sum_array = matrix(nrow = ncol(dat), 
#                      ncol = length(gene_list), 
#                      dimnames = list(colnames(dat),
#                                      names(gene_list)))
#   
#   lapply(names(gene_list), function(idx){
#     
#     sum_array[,idx] <<- apply(dat[rownames(dat) %in% gene_list[[idx]],], 2, sum)
#     
#   })
#   # return the sum array
#   return(sum_array)
#   
# }  


#### cell type enrichment ####

# function to assess the relative enrichment of cell type fractions compared to background
# specify matrix of cell type fractions and an index which labels each spot
# the background label name must also be specified
# mat.ct must have rows as spots and columns as features

cell_type_enrichment = function(mat.ct, spot_labels, background_label = 'other', enr_mode = 'empirical',
                                num_bg_samples = 1000, num_sample = 100){
  
  # define number of spots to sample
  num_niche = table(spot_labels[spot_labels != background_label])
  
  # for background I will use all spots

  sample_idx = replicate(num_bg_samples, sample(rownames(mat.ct), num_sample))
  
  niches = unique(spot_labels[spot_labels != background_label])
  if(enr_mode == 'empirical'){
    
    # calculate background frequencies
    # output has rows as cell types and random samples as columns
    bg_frac = apply(sample_idx, 2, function(spots){
      apply(mat.ct[spots,],2, mean)
    })
    bg_frac_mean = apply(bg_frac, 1, mean)
    niche_res = do.call(rbind, lapply(niches, function(niche){
      niche_mean = colMeans(mat.ct[spot_labels == niche,])
      # compute fraction of sample means that are greater than the niche mean
      emp = sapply(names(niche_mean), function(celltype){ sum(bg_frac[celltype,] > niche_mean[celltype])}) %>% as.matrix()
      # compute the log2 fold change between the mean of samples and the niche mean
      lfc = log2(niche_mean) - log2(bg_frac_mean[names(niche_mean)])
      dat_return = data.frame(celltype = names(niche_mean), 
                              niche = niche,
                              niche_mean_frac = niche_mean,
                              num_gt = emp,
                              log2FC = lfc[names(niche_mean)])
      return(dat_return)
    })
    )
  
  }
  niche_res = niche_res %>%
    mutate(emp_p = num_gt / num_bg_samples) %>%
    group_by(niche) %>%
    mutate(padj = p.adjust(emp_p, method = 'BH')) %>%
    ungroup()

  return(niche_res)
}


#### extractPositiveSpots ####

# function to take a binarized matrix for presence of certain cell types and extract the corresponding cells from another
# matrix with these

extractPositiveSpots = function(df, pos_mat){
  pos_mat = pos_mat[rownames(df), colnames(df)]
  idx = rownames(pos_mat)
  df_pos = do.call(rbind, lapply(colnames(df), function(ct){
    data.frame(celltype = ct, 
               value = df[pos_mat[,ct],ct], 
               barcode = idx[pos_mat[,ct]])
  })
  )
  return(df_pos)
}

#### remove_barcode ####
# function to remove the barcode component of a rowname and keeping the first part after the _
# which is typically the section index

remove_barcode = function(bc){
  gsub("_[AGCT]*-[0-9]+$", "", bc)
}


#### compute.background.level.df ####
# function to run the background level computation using dataframes and returning key parameters instead
# of having to go through spacefold object

compute.background.level.df = function (df.list, posterior.cutoff = 0.7, theta.cutoffs.user = NULL, 
          Znk.cutoffs.user = NULL, seed = NULL, num_mix = 2) 
{
  if (is.null(seed)) 
    set.seed(seed)
  Znk <- df.list[['Znk']]
  theta <- df.list[['theta']]
  theta.cutoffs <- c()
  cat("fitting mixture models on theta... \n")
  cat("Current cell type: ")
  for (i in 1:ncol(theta)) {
    cat(colnames(theta)[i], " ")
    capture.output({
      fit = gammamixEM(theta[, i], k = num_mix, maxit = 10000, 
                       maxrestarts = 100, mom.start = F)
    })
    cls <- apply(fit$posterior, 1, which.max)
    posterior <- fit$posterior
    cls.mean <- fit$gamma.pars["alpha", ] * fit$gamma.pars["beta", 
    ]
    # if (cls.mean[1] > cls.mean[2]) 
    #   max.cls <- 1
    # else max.cls <- 2
    # take second highest distribution
    cls.order = rank(cls.mean)
    max.cls = which(cls.order == 2)
    high.cls.idx <- posterior[, max.cls] > posterior.cutoff
    while (length(high.cls.idx) == 0) {
      posterior.cutoff <- posterior.cutoff - 0.05
      high.cls.idx <- posterior[, max.cls] > posterior.cutoff
    }
    theta.cutoff.i <- min(theta[high.cls.idx, i])
    theta.cutoffs <- c(theta.cutoffs, theta.cutoff.i)
  }
  names(theta.cutoffs) <- colnames(theta)
  cat("\n")
  if (!is.null(theta.cutoffs.user)) {
    if (length(theta.cutoffs.user) == 1) {
      theta.cutoffs.user <- rep(theta.cutoffs.user, length(theta.cutoffs))
      names(theta.cutoffs.user) <- names(theta.cutoffs)
    }
    stopifnot(!is.null(names(theta.cutoffs.user)))
    auto.cutoff.selected <- theta.cutoffs[names(theta.cutoffs.user)]
    user.higher.idx <- names(auto.cutoff.selected)[auto.cutoff.selected < 
                                                     theta.cutoffs.user]
    if (length(user.higher.idx) > 0) 
      theta.cutoffs[user.higher.idx] <- theta.cutoffs.user[user.higher.idx]
  }
  
  Znk.cutoffs <- c()
  cat("fitting mixture models on Znk... \n")
  cat("Current cell type: ")
  for (i in 1:ncol(Znk)) {
    cat(colnames(Znk)[i], " ")
    capture.output({
      fit <- Mclust(Znk[, i], G = num_mix)
    })
    posterior <- fit$z
    cls.mean <- fit$parameters$mean
    # if (cls.mean[1] > cls.mean[2]) 
    #   max.cls <- 1
    # else max.cls <- 2
    cls.order = rank(cls.mean)
    max.cls = which(cls.order == 2)
    high.cls.idx <- posterior[, max.cls] > posterior.cutoff
    while (length(high.cls.idx) == 0) {
      posterior.cutoff <- posterior.cutoff - 0.05
      high.cls.idx <- posterior[, max.cls] > posterior.cutoff
    }
    Znk.cutoff.i <- min(Znk[high.cls.idx, i])
    Znk.cutoffs <- c(Znk.cutoffs, Znk.cutoff.i)
  }
  names(Znk.cutoffs) <- colnames(Znk)
  cat("\n")
  if (!is.null(Znk.cutoffs.user)) {
    if (length(Znk.cutoffs.user) == 1) {
      Znk.cutoffs.user <- rep(Znk.cutoffs.user, length(Znk.cutoffs))
      names(Znk.cutoffs.user) <- names(Znk.cutoffs)
    }
    stopifnot(!is.null(names(Znk.cutoffs.user)))
    auto.cutoff.selected <- Znk.cutoffs[names(Znk.cutoffs.user)]
    user.higher.idx <- names(auto.cutoff.selected)[auto.cutoff.selected < 
                                                     Znk.cutoffs.user]
    if (length(user.higher.idx) > 0) 
      Znk.cutoffs[user.higher.idx] <- Znk.cutoffs.user[user.higher.idx]
  }
  
  # return list of cutoffs for Znk and theta
  rl = list(znk_cut = Znk.cutoffs, 
            theta_cut = theta.cutoffs)
  return(rl)
}

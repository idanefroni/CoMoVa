library("seqinr")
library("data.tree")
library("Biostrings")

#######################################################################################
# Setup global variables
#######################################################################################

# The maximum number of candidate orthologs for each gene
max_orthologs = 8

# Species lists
# When adding a new genome. Add here.
#
# dicots, monocots, brasicacea

dicots = c("Athaliana","Acoerulea","Ahalleri","Ahypochondriacus" ,
  "Alyrata"           ,"Boleraceacapitata" ,"BrapaFPsc"         ,"Bstricta",
  "Cclementina"       ,"Cgrandiflora",         
  "Cpapaya"           ,"Crubella"          ,"Csativus"          ,"Csinensis",        
  "Dcarota"           ,"Egrandis"          ,"Esalsugineum"      ,"Fvesca",           
  "Gmax"              ,"Graimondii"        ,"Kfedtschenkoi"     ,"Klaxiflora",       
  "Lusitatissimum"    ,"Mdomestica"        ,"Mguttatus"         ,"Mtruncatula",      
  "Mesculenta"        ,"Ptrichocarpa"      ,"Pvulgaris"         ,"Ppersica",
  "Rcommunis",        "Slycopersicum"      ,"Spurpurea",        
  "Stuberosum"        ,"Tcacao"            ,"Tpratense"         ,"Vvinifera")
monocots = c("Acomosus","Bdistachyon","Bstacei","Osativa","Sbicolor","Sitalica","Sviridis","Zmays")
brasicacea = c("Athaliana","Ahalleri",
               "Alyrata"           ,"Boleraceacapitata" ,"BrapaFPsc",         
               "Bvulgaris"            ,"Cgrandiflora"      ,"Chirsuta",         
               "Crubella","Esalsugineum")


#####################################################################################################
# Building the ortholog tables from reciprocal blast tables
#####################################################################################################

####### remove everything after the last dot
#
stripDot <- function(x, dot="\\.", trim=1) {
	if( regexpr(dot,x)[[1]]>0) {
		substr(x,1, regexpr(dot,x)[[1]]-trim)
	 } else {
	    x
         }
}

######### eval2ranks - accept a table from blast (outfmt 6)
# 1. filter different parts of the same gene
# 2. remove more the max_ortholog orthologs (picking just the highest Eval ones)
# 3. Return a table of orthologs with ranks (1- best)
#

eval2ranks <- function(blastres) {
  genenames = unique(blastres[,1])
  blastres=blastres[,1:3]
  
  # for gene parts, pick the best one
  do.call(rbind, 
    lapply(genenames, function (g) {
      blastresg = blastres[blastres[,1]==g,]
    
      #pick just the best orth (always the first)
      blastresg = blastresg[!duplicated(blastresg[,2]),]
      # compute ranks
      blastresg[,3] = seq(1, nrow(blastresg),1)
      # filter max hits
      
      blastresg = blastresg[ blastresg[,3] <= max_orthologs,]
  }))
  
  
}

########## buildOrthologs - build an ortholog table from reciprocal blast
# a2b - eval2rank output table of blast species A vs. species B
# b2a - eval2rank output table of blast species B vs. species A
# gene names must match
#
# returns a table of orthologs with the first column is ortholog quality score (2-101) and list of up to max_orthologs orthologs
#     for each gene
#
#

buildOrthologs <- function(a2b, b2a) {
  # Relaxation criteria
  relax = 4
  
  genenames_a= unique(a2b[,1])
  ortho_mat = matrix(nrow=length(genenames_a), ncol= (max_orthologs+1), data="")
  rownames(ortho_mat)= genenames_a
  
  #make genematrix
  for(g in genenames_a) {
    putative_orth = (a2b[a2b[,1]==g,1:3])
    if(!is.matrix(putative_orth)) {
      putative_orth = matrix(nrow=1, ncol=3,data = putative_orth) 
    }
    putative_orth = cbind(putative_orth,t(unlist(apply(putative_orth, 1, 
              function(x) {
                # handle two cases. One, we have duplicates. Others, we don't.
                orths = b2a[b2a[,1]==x[2] & b2a[,2]==g,c(1,3)]
                if( is.matrix(orths) && length(orths)>0){
                  orths= as.vector(orths[1,])
                }
                if(length( orths)==2) {
                  orths
                } else {
                  c("",100)
                } }))))
    putative_orth = cbind(putative_orth, as.numeric(putative_orth[,3])+ as.numeric(as.vector(putative_orth[,5])))
    best_orth_score = min(as.numeric(putative_orth[, ncol(putative_orth)]))
    num_best_orth = min(max_orthologs, sum(as.numeric(putative_orth[,6]) <= best_orth_score + relax))
    ortho_mat[g,2:(num_best_orth+1)]= head(putative_orth[as.numeric(putative_orth[,6]) <=best_orth_score + relax,2],n = num_best_orth)
    ortho_mat[g,1] = best_orth_score
  }
 ortho_mat
}

#####################################################################################################
##### Access ortholog table
#####################################################################################################

###### getOrthologs
#
# returns the orthologs for a give gene from an ortholog table
#
# strict - return just the first ortholog 
# allorth - return all candidate ortholog
# default behaviour( strict=FALSE, allorth=FALSE) is to return a number of orthologs based on the quality score
#   (for high score, return just 1 gene. For lower scores, return several)
#

getOrthologs <- function(locus, orthologs, strict=FALSE, allorth=FALSE) {
  if(!(locus %in% rownames(orthologs))) {
    return(NA)
  }
  
  orth = orthologs[locus,]
  if(as.numeric(orth[1])<=3 && !allorth) {
    a=orth[2]
  } else {
    a = unique(orth[2:7])
    a= a[a!=""]
  }
  if(strict) { a[1] } else {a}
}

###### getArabidopsisOrtholog
#
# returns the Arabidopsis ortholog for the a given gene
#
# strict - return just the first ortholog 
#

getArabidopsisOrtholog <- function(locus, orthologs, strict=TRUE) {
  if(strict) {
    arao = (rownames(orthologs)[orthologs[,2]==locus])[1]
    return(arao)
  }  
}

########## getAllOrthologs
#
#  returns the orthologs of a given gene for all species
#   see getOrthologs for strict\allorth
#

getAllOrthologs <- function(locus, orthologs, strict, allorth, orthologlist=names(orthologs)) {
  a=unlist(lapply(orthologs[orthologlist], function(x) { getOrthologs(locus, x, strict, allorth) } ))
  a[!is.na(a)]
}

########## getAllOrthologsbyRE 
#
# same as getAllOrthologs, but takes into account the similarity between the RE of the different orthologs
#  for each gene and species, will select of the candidate ortholog, the one with the most similar RE composition
# 

getAllOrthologsByRE <- function(locus, orthologs, RE, orthologlist=names(orthologs)) {
  baseRE = sort(getREVariants("Athaliana", locus, RE)[,3])
  
  a=unlist(lapply(orthologlist,
          function(x) { 
              orth = getOrthologs(locus, orthologs[[x]], strict=FALSE, allorth = TRUE)
              if(length(orth)>1) {
                bestorth = orth[1]
                bestRE = getREDist(baseRE, getREVariants(x , orth[1], RE)[,3])
                for(i in 2:length(orth)) {
                  
                  if(getREDist(baseRE, getREVariants(x, orth[i], RE)[,3])> bestRE) {
                    bestorth = orth[i]
                    bestRE = getREDist(baseRE, sort(getREVariants(x , orth[i], RE)[,3]))
                  }
                }
                bestorth
              } else {
                orth
              }
            
            }
          ))
  names(a) = orthologlist
  a[!is.na(a)]
}

####### getREDist 
#
# utility function for getAllOrthologsbyRE.
# returns the distance between two RE vectors
#

getREDist <- function (a,b) {
  dist = pmatch(a,b)
  sum(!is.na(dist))
}


#######################################################################################
#### Response Element Database
#######################################################################################


getREVariants <- function(species, locus, redb) {
    redb[[species]][ redb[[species]][,1] == locus, ]
}

##################################################################################################
#### Processing RE tables
##################################################################################################

######### loadcontable
#
# loads all conservation data files (from getConservationData)
# if we have background (usually the flanking sequence), fit distribution and calculate p-val

loadcontable <- function(name, background=NULL, directory="") {
  
  con = list()
  con$angio = as.matrix(read.csv(paste(directory,"angio_",name,".csv", sep=""), as.is=TRUE, row.names=1)) 
  con$dicots = as.matrix(read.csv(paste(directory,"dicots_",name,".csv", sep=""), as.is=TRUE, row.names=1))
  con$monocots = as.matrix(read.csv(paste(directory,"monocots_",name,".csv", sep=""), as.is=TRUE, row.names=1))
  con$bra = as.matrix(read.csv(paste(directory,"bra_",name,".csv", sep=""), as.is=TRUE, row.names=1))
  con$rosid = as.matrix(read.csv(paste(directory,"rosid_",name,".csv", sep=""), as.is=TRUE, row.names=1))
  
  if(!is.null(background)) {
    con$angio_pval = makePvalMatrix(con$angio, background$angio)
    con$dicots_pval = makePvalMatrix(con$dicots, background$dicots)
    con$monocots_pval = makePvalMatrix(con$monocots, background$monocots)
    con$bra_pval = makePvalMatrix(con$bra, background$bra)
    con$rosid_pval = makePvalMatrix(con$rosid, background$rosid)
    
    } else {
    con$angio_pval = makePvalMatrix(con$angio, con$angio)
    con$dicots_pval = makePvalMatrix(con$dicots, con$dicots)
    con$monocots_pval = makePvalMatrix(con$monocots, con$monocots)
    con$bra_pval = makePvalMatrix(con$bra, con$bra)
    con$rosid_pval = makePvalMatrix(con$rosid, con$rosid)
  }
  con
}

########### makePvalMatrix
#
# utility function that accepts a matrix of conservation scores, fits a gamma distribution to the background and returns a matrix of conservation pvalues
#

makePvalMatrix <- function(val, background) {
  pval = val
  for(i in intersect(colnames(background),colnames(val))) {
#    backlist = background[background[,i]>=0,i]
    backlist = background[,i]
    backlist[backlist<0]=0
    
    #fit a negative binomial dist
    fitted_dist = tryCatch({ 
      fitdistrplus::fitdist(backlist,"nbinom", discrete = TRUE) 
    }, error=function(e) {
      NULL
    })
      
    if(!is.null(fitted_dist)) {
      pval[,i] = sapply(val[,i], function(x) {
        if(x>0) {
          -pnbinom(x, fitted_dist$estimate[1], mu=fitted_dist$estimate[2], log.p=TRUE, lower.tail = FALSE)
        } else {
          0
        }})
    } else {
      pval[,i] = rep(0,nrow(val))
    }
  }
  pval
} 
##################################################################################################
#### Setting up tree
##################################################################################################


getTemplateTree <- function() {
  leaves = list()
  ang = Node$new("Angiosperm")
  mon = ang$AddChild("Monocots")
  leaves[["Acomosus"]] = mon$AddChild("Acomosus")
  grass = mon$AddChild("Grass")
  g1 = grass$AddChild("G1")
  g2 = g1$AddChild("G2")
  leaves[["Bdistachyon"]] = g2$AddChild("Bdistachyon")
  leaves[["Bstacei"]] = g2$AddChild("Bstacei")
  leaves[["Osativa"]] = g1$AddChild("Osativa")
  g3 = grass$AddChild("G3")
  pani = g3$AddChild("panicoideae")
  p1 = pani$AddChild("p1")
  p11 = p1$AddChild("p11")
  leaves[["Sitalica"]]= p11$AddChild("Sitalica")
  leaves[["Sviridis"]]= p11$AddChild("Sviridis")
  p2 = pani$AddChild("p2")
  leaves[["Zmays"]] = p2$AddChild("Zmays")
  leaves[["Sbicolor"]] = p2$AddChild("Sbicolor")
  eudi = ang$AddChild("eudictos")
  leaves[["Acoerulea"]] = eudi$AddChild("Acoerulea")
  penta = eudi$AddChild("pentapetalae")
  upasterids = penta$AddChild("upasterid")
  leaves[["Ahypochondriacus"]] = upasterids$AddChild("Ahypochondriacus")
  asterids = upasterids$AddChild("Asterids")
  leaves[["Dcarota"]] = asterids$AddChild("Dcarota")
  solu1 = asterids$AddChild("solu1")
  leaves[["Mguttatus"]] = solu1$AddChild("Mguttatus")
  sol = solu1$AddChild("Solanacea")
  leaves[["Slycopersicum"]] = sol$AddChild("Slycopersicum")
  leaves[["Stuberosum"]] = sol$AddChild("Stuberosum")
  rosup = penta$AddChild("rosup")
  kk = rosup$AddChild("kk")
  leaves[["Klaxiflora"]]=kk$AddChild("Klaxiflora")
  leaves[["Kfedtschenkoi"]] = kk$AddChild("Kfedtschenkoi")
  rosid = rosup$AddChild("rosid")
  leaves[["Vvinifera"]] = rosid$AddChild("Vvinifera")
  m2 = rosid$AddChild("m2")
  leaves[["Egrandis"]] = m2$AddChild("Egrandis")
  m1 = m2$AddChild("m1")
  fab = m1$AddChild("fabidea")
  fab1 = fab$AddChild("fabacea")
  f1 = fab1$AddChild("f1")
  f2 = fab1$AddChild("f2")
  leaves[["Gmax"]]= f1$AddChild("Gmax")
  leaves[["Pvulgaris"]] = f1$AddChild("Pvulgaris")
  leaves[["Mtruncatula"]] = f2$AddChild("Mtruncatula")
  leaves[["Tpratense"]] = f2$AddChild("Tpratense")
  fx1 = fab$AddChild("f1")
  leaves[["Csativus"]] = fx1$AddChild("Csativus")
  fx2 = fx1$AddChild("f2")
  leaves[["Fvesca"]] = fx2$AddChild("Fvesca")
  fx3 = fx2$AddChild("f3")
  leaves[["Mdomestica"]] = fx3$AddChild("Mdomestica")
  leaves[["Ppersica"]] = fx3$AddChild("Ppersica")
  malvidae = m1$AddChild("Malvidae")
  malpig = malvidae$AddChild("malpig")
  mx1 = malpig$AddChild("mx1")
  leaves[["Rcommunis"]] = mx1$AddChild("Rcommunis")
  leaves[["Mesculenta"]] = mx1$AddChild("Mesculenta")
  mz2 = malpig$AddChild("m2")
  leaves[["Lusitatissimum"]] = mz2$AddChild("Lusitatissimum")
  mz3 = m2$AddChild("m3")
  leaves[["Ptrichocarpa"]] = mz3$AddChild("Ptrichocarpa")
  leaves[["Spurpurea"]] = mz3$AddChild("Spurpurea")
  sbm = malvidae$AddChild("SBM")
  citrus = sbm$AddChild("citrus")
  leaves[["Cclementina"]] = citrus$AddChild("Cclementina")
  leaves[["Csinensis"]] = citrus$AddChild("Csinensis")
  bramal = sbm$AddChild("Brassicales-Malvales")
  my1 = bramal$AddChild("m1")
  leaves[["Graimondii"]] = my1$AddChild("Graimondii")
  leaves[["Tcacao"]] = my1$AddChild("Tcacao")
  mq1 = bramal$AddChild("brassup")
  leaves[["Cpapaya"]] = mq1$AddChild("Cpapaya")
  brassica = mq1$AddChild("brassica")
  b1 = brassica$AddChild("b1")
  leaves[["Esalsugineum"]] = b1$AddChild("Esalsugineum")
  b2 = b1$AddChild("b2")
  leaves[["BrapaFPsc"]] = b2$AddChild("BrapaFPsc")
  leaves[["Boleraceacapitata"]] = b2$AddChild("Boleraceacapitata")
  bra1 = brassica$AddChild("bra1")
  leaves[["Bstricta"]] = bra1$AddChild("Bstricta")
  bra2 = bra1$AddChild("bra2")
  bra3 = bra2$AddChild("bra3")
  leaves[["Cgrandiflora"]] = bra3$AddChild("Cgrandiflora")
  leaves[["Crubella"]] = bra3$AddChild("Crubella")
  bra4= bra2$AddChild("bra4")
  leaves[["Athaliana"]] = bra4$AddChild("Athaliana")
  bra5 = bra4$AddChild("bra5")
  leaves[["Ahalleri"]] = bra5$AddChild("Ahalleri")
  leaves[["Alyrata"]] = bra5$AddChild("Alyrata")

  intnode = list(ang, mon,  grass,g1,g2, g3, pani, p1,p11, p2, eudi, penta, asterids, rosup, kk, rosid, m2,m1,fab,fab1,f1,f2, fx1,fx2,fx3,
                 malvidae, malpig,mx1, mz2, mz3, sbm, citrus, bramal,my1, mq1,b1,b2, bra1, bra2, bra3, bra4,bra5, solu1, sol)
  names(intnode)[[35]] = "Brassicaceae"
  names(intnode)[[31]]="sbm"
  names(intnode)[[1]] = "Angiosperms"
  names(intnode)[[2]] = "Monocots"
  names(intnode)[[3]] = "Grasses"
  names(intnode)[[11]] = "Dicots"
  names(intnode)[[16]] = "Rosids"
  viri = fillTree(ang, "")
  list(ang, leaves, intnode)
}

##################################################################################################
#### Processing RE conservation tree
##################################################################################################

##### loadTree
#
# build an ortholog variant tree for a give gene using the template tree "tree".
#

loadTree <- function(tree, gene, redb) {
  orth = c(gene, getAllOrthologsByRE(gene, orthologs, redb))
  names(orth)[1] = "Athaliana"

  for(i in intersect(names(orth), names(tree[[2]]))) {
    tree[[2]][[i]]$variant = unique(getREVariants( i, orth[i], redb)[,3])
  }

  tree[[1]] = populateTree(tree[[1]])
  tree = optimizeTree_Parsimony(tree,100)
  tree
}
##### duplicateTree
#
# returns a copy of the tree
#
duplicateTree <- function(source, target) {
   target[[1]] = int_copyTree(source[[1]], target[[1]])
   target
}

##### int_copyTree
#
# utility function used by duplicateTree
#
int_copyTree <- function(source, target) {
  target$variant = source$variant
  if(! (source$isLeaf)) {
    for(i in 1:length(source$children)) {
      target$children[[i]] = int_copyTree(source$children[[i]], target$children[[i]])
    }
  }
  return(target)
}

##### fillTree
#
# utility function to set the variant of all nodes in "tree" to "value"
#

fillTree <- function(tree, value) {
  tree$variant = value
  if(tree$isLeaf) {
    return(tree)
  } else {
    for(i in 1:length(tree$children)) {
      tree$children[[i]] = fillTree(tree$children[[i]], value)
    }
    return(tree)
  }
}

##### addToTree
#
# utility function to add a variant "varToAdd" to all internal nodes of "tree"
#

addToTree <- function(tree, varToAdd) {
  if(!tree$isLeaf) {
    for(i in 1:length(tree$children)) {
      tree$children[[i]] = addToTree(tree$children[[i]], varToAdd)
    }
    # remove inconsistent variants
    tree$variant = unique(c(tree$variant, varToAdd))
  }
  return(tree)
}

##### removeFromTree
#
# utility function to removes a variant "varToRemove" from all internal nodes of "tree", unless all children have the variant
#

removeFromTree <- function(tree, varToRemove) {
  if(!tree$isLeaf) {
    for(i in 1:length(tree$children)) {
      tree$children[[i]] = removeFromTree(tree$children[[i]], varToRemove)
    }
    # add inconsistent variants
    tree$variant = setdiff(tree$variant, varToRemove)
  }
  return(tree)
}

######## changes
#
# utility function to count the number of changes between two set of variants. Used by countParismony
#
changes <- function(a,b, variant) {
  if(variant=="") {
    length(union(a,b))- length(intersect(a,b))
  } else {  
    length( intersect( union(a,b), variant)) - length( intersect(variant, intersect(a,b)))
  }
}

######## countParsimony
#
# Compute the number of changes between variants occuring along a tree
#
countParsimony <- function(tree, variant="") {
  if(tree$isLeaf) {
    return (0)
  } else {
    c1 = changes(tree$variant, tree$children[[1]]$variant, variant ) + countParsimony(tree$children[[1]], variant)
    if(length(tree$children)>1) {
      c2 = changes(tree$variant, tree$children[[2]]$variant, variant ) + countParsimony(tree$children[[2]], variant)
    } else {
      c2=0
    }
  }
  c1+c2
}

######### trimTree
#
# scans the tree for internal nodes that are inconsistent with children and remove variant
#

trimTree <- function(tree, tight=FALSE, consultSib = FALSE) {
  if(!tree$isLeaf) {
    tree$children[[1]] = trimTree(tree$children[[1]], tight, consultSib)
    if(length(tree$children)>1) {
      tree$children[[2]] = trimTree(tree$children[[2]], tight, consultSib)
      if(tight) {
        childvariants = intersect(tree$children[[1]]$variant, tree$children[[2]]$variant )
      } else if (consultSib && !is.null(tree$parent)) {
          sib1 = tree$parent$children[[1]]
          sib2 = tree$parent$children[[1]]
          childvariants = c(tree$children[[1]]$variant, tree$children[[2]]$variant)
          childvariants = unique(c(intersect(childvariants, sib1$variant),
                            intersect(childvariants, sib2$variant),
                            intersect(tree$children[[1]]$variant, tree$children[[2]]$variant)))
          childvariantshelps = rep( FALSE ,length(childvariants))
          tree$variant = childvariants

          if(length(childvariants)>0) {
            # now check all the possible parsimonies and see if deleting helps
            for(i in 1:length(childvariants)) {
              cur_parsimony = countParsimony(tree, childvariants[i])
              tree$variant = setdiff(tree$variant, childvariants[i])
              new_parsimony = countParsimony(tree, childvariants[i])
              if(new_parsimony >= cur_parsimony) {
                childvariantshelps[i] = TRUE
              }
            }
            childvariants = childvariants[childvariantshelps]
          }
      } else {
        childvariants = c(tree$children[[1]]$variant, tree$children[[2]]$variant )
      }
    } else {
      childvariants = tree$children[[1]]$variant
    }
    tree$variant = childvariants
  }
  return(tree)
}

######### plasterTree
#
# scans the tree for internal nodes that are inconsistent with children and add variant
#

plasterTree <- function(tree, tight=FALSE) {
  if(!tree$isLeaf) {
    childvariants = c()
    for(i in 1:length(tree$children)) {
      tree$children[[i]] = plasterTree(tree$children[[i]],tight)
      childvariants = c(childvariants, tree$children[[i]]$variant)
    }
    # add inconsistent variants
    tree$variant = unique(c(tree$variant, unique(childvariants[duplicated(childvariants)])))
  }
  return(tree)
}

######### populateTree
#
# fill the tree internal nodes based on the children's variants
#

populateTree <- function(tree, hangon=c("")) {
  if(tree$isLeaf) {
    return (tree)
  } else {
    c1 = populateTree(tree$children[[1]])
    tree$variant = c1$variant

    if(length(tree$children)>1) {
        c2 = populateTree(tree$children[[2]])
        priority = intersect(hangon, c(c1$variant, c2$variant))
        jointvar = setdiff( union(priority,intersect(tree$variant, c2$variant)),c(""))
        if(length(jointvar)==0) {
          jointvar = setdiff(union(tree$variant, c2$variant),c(""))
        }
        tree$variant = jointvar
    }

    return(tree)
  }
}

######### optimizeTree_Parsimony
#
# find and returns a maximum parsimony tree for a given set of loaded leaves.
#  it's a bit of a hack and there are better ways to do it, but it works.
#

optimizeTree_Parsimony <- function(tree, depth=3, numofoptimizationsteps=2) {
  vartable = sort(table(unlist(sapply(tree[[2]], function(x) { x$variant}))), decreasing = TRUE)
  
  maxParismony = countParsimony(tree[[1]])

  depth=min(depth, length(vartable)/2)
  
  #optimize tree by adding and trimming
 
  for(n in 1:numofoptimizationsteps) {
    for(i in 1:depth) {
      temp_tree = getTemplateTree()
      temp_tree = duplicateTree(tree, temp_tree)
      if(n==1 || n==numofoptimizationsteps) {
        nodeToOptimize = temp_tree[[1]]
      } else {
        nodeToOptimize = sample(temp_tree[[3]],1)
      }
      #remove and plaster
      nodeToOptimize = removeFromTree(nodeToOptimize, names(vartable)[length(vartable)-i])
      nodeToOptimize = plasterTree(nodeToOptimize)
      temp_parsimony = countParsimony(temp_tree[[1]])
      if(temp_parsimony< maxParismony) {
        tree = duplicateTree(temp_tree, tree)
        maxParismony = temp_parsimony
      } else {
        temp_tree = duplicateTree(tree, temp_tree)
      }
    
      #add and trim
      nodeToOptimize = addToTree(nodeToOptimize, names(vartable)[i])
      nodeToOptimize = trimTree(nodeToOptimize,consultSib = TRUE)
      temp_parsimony = countParsimony(temp_tree[[1]])
      if(temp_parsimony< maxParismony) {
        tree = duplicateTree(temp_tree, tree)
        maxParismony = temp_parsimony
      } else {
        temp_tree = duplicateTree(tree, temp_tree)
      }

      #tighten
      nodeToOptimize = trimTree(nodeToOptimize, consultSib = TRUE)
      temp_parsimony = countParsimony(temp_tree[[1]])
      
      if(temp_parsimony< maxParismony) {
        tree = duplicateTree(temp_tree, tree)
        maxParismony = temp_parsimony
      } else {
        temp_tree = duplicateTree(tree, temp_tree)
      }
    }
    temp_tree = getTemplateTree()
    temp_tree = duplicateTree(tree, temp_tree)
    temp_tree[[1]] = trimTree(temp_tree[[1]], tight=TRUE)
    if(countParsimony(temp_tree[[1]]) < countParsimony(tree[[1]])) {
      tree = duplicateTree(temp_tree, tree)
    }
  }
  tree[[1]]$variant = unique(tree[[1]]$variant)
  return(tree)
}

##################################################################################################
#### Calculate conservation score
##################################################################################################

########## countVariantInTree
#
# return how many leaves have the variant
#

countVariantInTree <- function(tree, variant) {
  if(tree$isLeaf) {
    if(variant %in% tree$variant ) { return(1)}  else { return(0)}
  }
  
  varchil = 0
  for(i in 1:length(tree$children)) { varchil = varchil + countTree(tree$children[[i]], variant)}
  
  return(varchil)
}

########## countTree
#
# utility function for countVariantInTree
#

countTree <- function(tree, variant) {
  if(tree$isLeaf) {
    if(variant %in% tree$variant) { return(1)} else { return(0)}
  } else {
    
    varchil=0
    for(i in 1:length(tree$children)) { varchil = varchil+ countTree(tree$children[[i]], variant)}
    return(varchil)
  }
}

########## getConservationScore
#
# return the conservation score (how many leaves have the variants - how many changes along the tree)
#

getConservationScore <- function(tree, variant, subTree="Angiosperms") {
  unlist(sapply(variant, function (x) { countVariantInTree(tree[[3]][[subTree]],x) - countParsimony(tree[[3]][[subTree]],x) }))
}

########## getParsimonyScore
#
# return the number of changes along the tree
#

getParismonyScore <- function(tree, variant, subTree="Angiosperms") {
  unlist(sapply(variant, function (x) { countParsimony(tree[[3]][[subTree]],x) }))
}
########## getCountScore
#
# return the number of variants along the tree
#

getCountScore <- function(tree, variant, subTree="Angiosperms") {
  unlist(sapply(variant, function (x) { countVariantInTree(tree[[3]][[subTree]],x) }))
}

########## getConservationData
#
# For a given gene and RE database - construct the tree and return the conservation score for each variant as a list with a slot for each clades:
#	angiosperms, monocots, dicots, brassica, rosids
#

getConservationData <- function(gene,redb) {
  if(!exists("variantDic")) {
    variantDic = sort(unique(redb[[1]][,3]))
  }
  
  t1 <- getTemplateTree()
  t1 <- loadTree(t1, gene, redb = redb)
  
  list(name = gene, 
       conserv = getConservationScore(t1, variantDic),
       monocots_conserv = getConservationScore(t1, variantDic, "Monocots"),
       dicots_conserv = getConservationScore(t1, variantDic,"Dicots"),
       brassica_conserv = getConservationScore(t1, variantDic, "Brassicaceae"),
       rosid_conserv = getConservationScore(t1, variantDic, "Rosids"))
}


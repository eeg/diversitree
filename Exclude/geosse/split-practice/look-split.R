### 27 Feb 2012
### Dan Rabosky asked how the E's work when different clades have different
### sampling.f.  Here, I'm working through how diversitree handles this.  Key
### answer: at a breakpoint between partitions, E's are switched so that they 
### reflect the parameters that will be used further down the tree.


library(diversitree)

tree <- read.tree(text="(tipA:1.5, tipB:1.5)nodeG;")
states <- c("tipA"=0, "tipB"=1)
lnL <- make.bisse(tree, states, sampling.f=0.9)
pars <- c(10, 20, 2, 3, 0.2, 0.2)
lnL(pars, intermediates=TRUE)
# attr(,"intermediates")$init
#          [,1] [,2]       [,3]
#     [1,]  0.1  0.1 0.19879571
#     [2,]  0.1  0.1 0.15056776
#     [3,]  0.9  0.0 9.54222855
#     [4,]  0.0  0.9 0.01068918
# attr(,"intermediates")$base
#                [,1]       [,2] [,3]
#     [1,] 0.19879571 0.19879571   NA
#     [2,] 0.15056776 0.15056776   NA
#     [3,] 0.97816026 0.97552814   NA
#     [4,] 0.02183974 0.02447186   NA
# attr(,"intermediates")$lq
#     [1] -12.64143 -16.40364   0.00000
# attr(,"intermediates")$root.p
#     [1] 0.998881056 0.001118944
# attr(,"vals")
#     [1] 0.19879571 0.15056776 9.54222855 0.01068918

tree <- read.tree(text="(tipC:2.25, tipD:2.25)nodeF;")
states <- c("tipC"=0, "tipD"=1)
lnL <- make.bisse(tree, states, sampling.f=0.4)
pars <- c(10, 20, 2, 3, 0.2, 0.2)
lnL(pars, intermediates=TRUE)
# attr(,"intermediates")$init
#          [,1] [,2]       [,3]
#     [1,]  0.6  0.6 0.19879611
#     [2,]  0.6  0.6 0.15056777
#     [3,]  0.4  0.0 9.43911489
#     [4,]  0.0  0.4 0.01521846
# attr(,"intermediates")$base
#                [,1]       [,2] [,3]
#     [1,] 0.19879611 0.19879611   NA
#     [2,] 0.15056777 0.15056777   NA
#     [3,] 0.97843395 0.96471662   NA
#     [4,] 0.02156605 0.03528338   NA
# attr(,"intermediates")$lq
#     [1] -18.03339 -21.40380   0.00000
# attr(,"intermediates")$root.p
#     [1] 0.998390319 0.001609681
# attr(,"vals")
#     [1] 0.19879611 0.15056777 9.43911489 0.01521846

tree <- read.tree(text="((tipA:1.5, tipB:1.5)nodeG:2.0, (tipC:2.25, tipD:2.25)nodeF:1.25)nodeE;")
states <- c("tipA"=0, "tipB"=1, "tipC"=0, "tipD"=1)
lnL <- make.bisse.split(tree, states, nodes=c("nodeG"), split.t=c(2.5), sampling.f=list(0.9, 0.4))

pars <- rep(pars, 2)
names(pars) <- argnames(lnL)

fixInNamespace("all.branches.split", "diversitree")
lnL(pars)


### C+D, at nodeF:
# Browse[2]> obj
# $init
#          [,1] [,2]        [,3]
#     [1,]  0.6  0.6 0.198799719
#     [2,]  0.6  0.6 0.150567852
#     [3,]  0.4  0.0 9.580926387
#     [4,]  0.0  0.4 0.008961169
# $base
#                [,1]       [,2] [,3]
#     [1,] 0.19879972 0.19879972   NA
#     [2,] 0.15056785 0.15056785   NA
#     [3,] 0.97816127 0.97948331   NA
#     [4,] 0.02183873 0.02051669   NA
# $lq
#     [1] -11.85642 -15.31385   0.00000

### C+D, but only go as far as nodeG
# Browse[2]> res[[2]]
# $lq
#     [1] -33.1085
# $base
#                [,1]
#     [1,] 0.19879610
#     [2,] 0.15056777
#     [3,] 0.97816303
#     [4,] 0.02183697
# $intermediates = obj

### C+D at nodeG's time (j=2 is daughter clade now; i = 1)
# Browse[2]> yj
#                [,1]
#     [1,] 0.19879610
#     [2,] 0.15056777
#     [3,] 0.97816303
#     [4,] 0.02183697
### run to breakpoint time = 2.5, using i's pars (but still j's samp.f)
# Browse[2]> tmp
#     [[1]]  # lq
#         [1] -8.219612
#     [[2]]  # base
#                   [,1]
#         [1,] 0.1987961
#         [2,] 0.1505678
#         [3,] 0.9781604
#         [4,] 0.0218396
### preset is the same as tmp, but also has $target = 1

### using A+B pars and also A+B's samp.f
### A+B, then down to nodeE from both sides
# Browse[2]> obj
#     $init
#              [,1] [,2] [,3]       [,4]      [,5]
#         [1,]   NA  0.1  0.1 0.19879610 0.1987961
#         [2,]   NA  0.1  0.1 0.15056777 0.1505678
#         [3,]   NA  0.9  0.0 9.56789981 0.5629213
#         [4,]   NA  0.0  0.9 0.00954284 0.4166691
#     $base
#                   [,1]       [,2]       [,3] [,4]       [,5]
#         [1,] 0.1987961 0.19879610 0.19879610   NA 0.19879610
#         [2,] 0.1505678 0.15056777 0.15056777   NA 0.15056777
#         [3,] 0.9781604 0.97789402 0.05756465   NA 0.97815244
#         [4,] 0.0218396 0.02210598 0.94243535   NA 0.02184756
#     $lq
#         [1]  -8.219612 -18.796561 -21.439317   0.000000 -10.811095

# Browse[2]> res[[1]]
#     $lq
#         [1] -59.26659
#     $base
#         [1] 0.19879610 0.15056777 9.56789981 0.00954284
#     $intermediates$init
#              [,1] [,2] [,3]       [,4]      [,5]
#         [1,]   NA  0.1  0.1 0.19879610 0.1987961
#         [2,]   NA  0.1  0.1 0.15056777 0.1505678
#         [3,]   NA  0.9  0.0 9.56789981 0.5629213
#         [4,]   NA  0.0  0.9 0.00954284 0.4166691
# 
#     $intermediates$base
#                   [,1]       [,2]       [,3] [,4]       [,5]
#         [1,] 0.1987961 0.19879610 0.19879610   NA 0.19879610
#         [2,] 0.1505678 0.15056777 0.15056777   NA 0.15056777
#         [3,] 0.9781604 0.97789402 0.05756465   NA 0.97815244
#         [4,] 0.0218396 0.02210598 0.94243535   NA 0.02184756
# 
#     $intermediates$lq
#         [1]  -8.219612 -18.796561 -21.439317   0.000000 -10.811095



tree <- read.tree(text = "((tip1:1,tip2:1)n5:1,tip3:2)n4;")
states <- c("tip1"=0, "tip2"=0, "tip3"=1)

lnL0 <- make.bisse(tree, states, sampling.f=c(0.1,0.9))
lnL <- make.bisse.split(tree, states, nodes=c("n5"), split.t=c(0), sampling.f=list(c(0.1,0.9), c(0.7,0.3)))

pars0 <- c(10, 20, 0.5, 0.5, 0.5, 0.5)
pars <- c(pars0, pars0)
names(pars0) <- argnames(lnL0)
names(pars) <- argnames(lnL)

lnL0(pars0)
lnL(pars)


fixInNamespace("all.branches.split", "diversitree")
lnL <- make.bisse.split(tree, states, nodes=c("n5"), split.t=c(0), sampling.f=list(c(0.1,0.9), c(0.7,0.3)))
lnL(pars)


Browse[1]> res
[[1]]
# still computing this
NULL

[[2]]
[[2]]$lq
[1] -17.11201

[[2]]$base
# before the node join, with logcomp removed from D
            [,1]
[1,] 0.048800673    # E0
[2,] 0.025595863    # E1
[3,] 0.995003492    # D0
[4,] 0.004996508    # D1

[[2]]$intermediates
[[2]]$intermediates$init
# after the node join
     [,1] [,2]       [,3]
[1,]  0.3  0.3 0.04880067
[2,]  0.7  0.7 0.02559586
[3,]  0.7  0.7 9.06843129
[4,]  0.0  0.0 0.04553802

[[2]]$intermediates$base
# before the node join
#    branch 1   branch 2
           [,1]       [,2] [,3]
[1,] 0.04880067 0.04880067   NA
[2,] 0.02559586 0.02559586   NA
[3,] 0.95228311 0.95228311   NA
[4,] 0.04771689 0.04771689   NA

[[2]]$intermediates$lq
[1] -9.66091 -9.66091  0.00000

# in all.branches.split
yj : for each daughter partition
    E0, E1 from parent partition, after time = this/daughter partition length
    D0, D1 for this partition
tmp : carry on down the rest of the branch (to next daughter depth or daughter-parent join)





$cache[[1]]
-----------
$cache[[1]]$tip.label
    [1] "tip3"
$cache[[1]]$node.label
    [1] "n4"
$cache[[1]]$len
    [1]  1  2 NA
$cache[[1]]$children
         [,1] [,2]
    [1,]   NA   NA
    [2,]   NA   NA
    [3,]    1    2
$cache[[1]]$parent
    [1]  3  3 NA
$cache[[1]]$order
    [1] 3
$cache[[1]]$root
    [1] 3
$cache[[1]]$n.tip
    [1] 1
$cache[[1]]$n.node
    [1] 1
$cache[[1]]$tips
    [1] 2
$cache[[1]]$height
    tip1 tip3   n4 
       1    2    0 
$cache[[1]]$depth
    tip1 tip3   n4 
       1    0    2 
$cache[[1]]$ancestors
    $cache[[1]]$ancestors[[1]]
        [1] 3
    $cache[[1]]$ancestors[[2]]
        [1] 3
    $cache[[1]]$ancestors[[3]]
        NULL
$cache[[1]]$edge
         [,1] [,2]
    [1,]    3    1
    [2,]    3    2
$cache[[1]]$edge.length
    [1] 1 2
$cache[[1]]$trailing.t0
    [1] 2
$cache[[1]]$daughters
    tip1 
       2 
$cache[[1]]$daughters.i
    [1] 1
$cache[[1]]$ny
    [1] 4
$cache[[1]]$k
    [1] 2
$cache[[1]]$tip.state
    tip3 
       1 
$cache[[1]]$sampling.f
    [1] 0.1 0.9
$cache[[1]]$y
    $cache[[1]]$y[[1]]
        $cache[[1]]$y[[1]]$y
            [1] 0.9 0.1 0.0 0.9
$cache[[1]]$y[[1]]$y.i
    [1] 2
$cache[[1]]$y[[1]]$target
    [1] 2
$cache[[1]]$y[[1]]$t.uniq
    [1] 2
$cache[[1]]$y[[1]]$unpack
    [1] 1

$cache[[2]]
-----------
$cache[[2]]$tip.label
    [1] "tip1" "tip2"
$cache[[2]]$node.label
    [1] "n5"
$cache[[2]]$len
    [1]  1  1 NA
$cache[[2]]$children
         [,1] [,2]
    [1,]   NA   NA
    [2,]   NA   NA
    [3,]    1    2
$cache[[2]]$parent
    [1]  3  3 NA
$cache[[2]]$order
    [1] 3
$cache[[2]]$root
    [1] 3
$cache[[2]]$n.tip
    [1] 2
$cache[[2]]$n.node
    [1] 1
$cache[[2]]$tips
    [1] 1 2
$cache[[2]]$height
    tip1 tip2   n5 
       1    1    0 
$cache[[2]]$depth
    tip1 tip2   n5 
       0    0    1 
$cache[[2]]$ancestors
    $cache[[2]]$ancestors[[1]]
        [1] 3
    $cache[[2]]$ancestors[[2]]
        [1] 3
    $cache[[2]]$ancestors[[3]]
        NULL
$cache[[2]]$edge
         [,1] [,2]
    [1,]    3    1
    [2,]    3    2
$cache[[2]]$edge.length
    [1] 1 1
$cache[[2]]$trailing.len
    n5 
     0 
$cache[[2]]$trailing.t0
    [1] 1
$cache[[2]]$ny
    [1] 4
$cache[[2]]$k
    [1] 2
$cache[[2]]$tip.state
    tip1 tip2 
       0    0 
$cache[[2]]$sampling.f
    [1] 0.5 0.5
$cache[[2]]$y
    $cache[[2]]$y[[1]]
        $cache[[2]]$y[[1]]$y
            [1] 0.3 0.7 0.7 0.0
$cache[[2]]$y[[1]]$y.i
    [1] 1
$cache[[2]]$y[[1]]$target
    [1] 1 2
$cache[[2]]$y[[1]]$t.uniq
    [1] 1
$cache[[2]]$y[[1]]$unpack
    [1] 1 1

overall
-------
$tip.tr
    tip1 tip2 tip3 
       2    2    1 
$parents
    [1] NA  1
$daughters
    $daughters[[1]]
        tip1 
           2 
    $daughters[[2]]
        NULL
$n.parts
    [1] 2
$order.parts
    [1] 2 1
$desc.parts
    $desc.parts[[1]]
        tip1 
           2 
    $desc.parts[[2]]
        NULL
$sampling.f
    $sampling.f[[1]]
        [1] 0.1 0.9
    $sampling.f[[2]]
        [1] 0.5 0.5
$aux.i
    [1] 1 2

# 
# edgeL = NULL
# edgeL = dlply(edges, "fromNode", function(edge){
#   #print(edge$fromNode[1])
#   a = NULL
#   a$edges = edge$toNode
#   a$weights = edge$weight
#   #print(a)
#   return(a)
#   # print(d) [[edge$fromNode[1]]]
# })
# 
# edgeTo = dlply(edges, "toNode", function(edge){
#   a = NULL
#   a$edges = edge$toNode
#   a$weights = edge$weight
#   return(a)
# })
# 
# g1 = new("graphNEL", nodes=as.character(nodes$id), edgeL=edgeL)
# 
# 
# layout = layout_with_fr(g, niter=15)
# plot.igraph(g, layout=layout, vertex.label=NA, vertex.size=1, 
#             vertex.color=nodes$color, edge.arrow.size=0)
# 
# layout = layout_with_lgl(g, maxiter=500)
# plot.igraph(g, layout=layout, vertex.label=NA, vertex.size=1, 
#             vertex.color=nodes$color, edge.arrow.size=0)
# 
# layout = layout_with_graphopt(g)
# plot.igraph(g, layout=layout, vertex.label=NA, vertex.size=1, 
#             vertex.color=nodes$color, edge.arrow.size=0)
#  
# 
# layout = layout_nicely(g)
# plot.igraph(g, layout=layout, vertex.label=NA, vertex.size=1, 
#             vertex.color=nodes$color, edge.arrow.size=0)
# 
# newEdges <- edges %>% 
#   rowwise() %>%
#   mutate(
#        fromModule = getModule(nodes, fromNode),
#        toModule = getModule(nodes, toNode)
#   )
# 
# getModule <- function(nodes, geneName){
#   return(filter(nodes, id == geneName)$color)
# }
# 
# filter(nodes, id == "AASDH"  )$color
# 
# a = ddply(edges, )
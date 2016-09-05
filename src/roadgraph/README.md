Class: MapGraph

Modifications made to MapGraph (what and why):

a) Introduced a map to maintain the mapping between the geographical
   location of an intersection, to our internal mapping of the vertex
   (MapVertex class), as described below. This helps in 0(1) lookup of
   the MapVertex class, given its geographical location.

b) We also introduced counters for tracking number of vertices and edges,
   so that we can look them up in 0(1) time complexity.

c) Introduced a validateEdge() private method which is invoked before we
   add an edge to the graph. It ensures that don't add a dangling edge in
   the graph, and also the edge properties like length, roadName and roadTypes
   are well defined.

------------------------------------------------------------

Class name: MapVertex

This class was introduced to encapsulate all the necessary properties of
an intersection. This includes its geographical location, a list of
outbound road segments, and also a reference to the parent
intersection, if this intersection is reached while doing a breath first
search of the city map. This structure ensures that we are able to reach all
the neighboring intersections without traversing the whole city. 

------------------------------------------------------------

Class name: MapEdge

This class stores all the attributes of the road segment between 2
intersections. The breath first algorithm of the city map works based upon
these attributes. This class is mainly a POJO model, with getter/setters
of the necessary attributes.

-------------------------------------------------------------

Overall Design Justification (4-6 sentences):

The main logic of this class sits in the bfs() function. The algorithm
followed for doing the BFS is very similar to what we discussed in the
lecture. Using the neighbor list rooted at each vertex, we travel the
graph breath wise, till we reach the destination intersection. The breath wise
traversal of the graph also ensures that the distance between source and
the destination is the shortest one.
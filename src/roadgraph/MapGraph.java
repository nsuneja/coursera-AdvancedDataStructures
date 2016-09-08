/**
 * @author UCSD MOOC development team and YOU
 * This class represents a graph of geographic locations.
 * Nodes in the graph are the intersections, and the edges are the road
 * segments.
 * 
 */
package roadgraph;


import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.function.Consumer;

import geography.GeographicPoint;
import util.GraphLoader;


/**
 * This class defines an edge in the graph.
 * @author nishantsuneja
 *
 */
class MapEdge
{
    private GeographicPoint from;
    private GeographicPoint to;
    private String roadName;
    private String roadType;
    private double length;

    public MapEdge(GeographicPoint from,
                   GeographicPoint to,
                   String roadName,
                   String roadType,
                   double length)
    {
        this.from = from;
        this.to = to;
        this.roadName = roadName;
        this.roadType = roadType;
        this.length = length;
    }

    public GeographicPoint getFrom()
    {
        return from;
    }
    public void setFrom(GeographicPoint from)
    {
        this.from = from;
    }
    public GeographicPoint getTo()
    {
        return to;
    }
    public void setTo(GeographicPoint to)
    {
        this.to = to;
    }
    public String getRoadName()
    {
        return roadName;
    }
    public void setRoadName(String roadName)
    {
        this.roadName = roadName;
    }
    public String getRoadType()
    {
        return roadType;
    }
    public void setRoadType(String roadType)
    {
        this.roadType = roadType;
    }
    public double getLength()
    {
        return length;
    }
    public void setLength(double length)
    {
        this.length = length;
    }

    @Override
    public String toString()
    {
        return "MapEdge [start=" + from + ", end=" + to + ", roadName="
                + roadName + ", roadType=" + roadType + ", length=" + length
                + "]";
    }

}


/**
 * This class defines a vertex in the graph. The class maintains the
 * following:
 * a)  List of all the outgoing edges from this vertex.
 * b)  The geographical location of the vertex.
 * c)  The geographical location of the parent of the vertex discovered as
 *     we traverse the graph.
 * 
 * @author nishantsuneja
 *
 */
class MapVertex
{
    private final GeographicPoint location;
    // This map stores the parent location of this vertex, while
    // traveling to a given destination vertex. The destination vertex
    // is the key and the value is the parent location of
    // this vertex, with respect to the given destination.
    private Map<GeographicPoint, GeographicPoint> parentLocationMap;
    private List<MapEdge> outEdges;
    private double distanceFromOrigin;

    public MapVertex(GeographicPoint location) {
        this.location = location;
        this.parentLocationMap =
                new HashMap<GeographicPoint, GeographicPoint>();
        this.outEdges = new ArrayList<MapEdge>();
        this.distanceFromOrigin = Double.MAX_VALUE;
    }

    public GeographicPoint getLocation()
    {
        return location;
    }
    public List<MapEdge> getOutEdges()
    {
        return outEdges;
    }
    public void setOutEdges(List<MapEdge> outEdges)
    {
        this.outEdges = outEdges;
    }

    /**
     * This vertex can have multiple parents based upon the end location
     * of the search. This method retrieves the parent location for the
     * given end location.
     * @param endLocation
     * @return
     */
    public GeographicPoint getParentLocation(GeographicPoint endLocation)
    {
        // The parent location of a vertex is now with respect to the
        // end vertex location.
        // NOTE: We use this method only after finding a path from source
        // to destination. At this point each node's (in the path) parent
        // should have been set.
        assert(parentLocationMap.containsKey(endLocation));
        return parentLocationMap.get(endLocation);
    }

    /**
     * This vertex can have multiple parents based upon the end location
     * of the search. This method sets the parent location for the
     * given end location.
     * @param parentLocation
     * @param endLocation
     */
    public void setParentLocation(GeographicPoint parentLocation,
                                  GeographicPoint endLocation)
    {
        assert((parentLocation != null) && (endLocation != null));
        parentLocationMap.put(endLocation, parentLocation);
    }

    public double getDistanceFromOrigin()
    {
        return distanceFromOrigin;
    }

    public void setDistanceFromOrigin(double distanceFromOrigin)
    {
        this.distanceFromOrigin = distanceFromOrigin;
    }

    /**
     * This method returns the edge length between this vertex and 
     * the given neighbor.
     * @param vertex
     * @return
     */
    public double getNeighborEdgeLength(MapVertex vertex) {
        for (MapEdge edge: outEdges) {
            if (vertex.getLocation().equals(edge.getTo())) {
                // Found the vertex
                return edge.getLength();
            }
        }

        // Unable to find the given vertex in the neighbor list
        return Double.MAX_VALUE;
    }

    /**
     * This method adds a neighbor to a vertex. If the same neighbor
     * already exists, we return false. Otherwise, true.
     * @param edge
     * @return
     */
    public boolean addOutEdge(MapEdge edge) {
        if (edge == null || outEdges.contains(edge)) {
            // The edge already exists.
            return false;
        }

        // Add the edge to the list
        outEdges.add(edge);

        return true;
    }

    @Override
    public String toString()
    {
        return "MapVertex [location=" + location + "]";
    }

    @Override
    public int hashCode()
    {
        return new Double(41*location.getX() + 43*location.getY()).intValue();
    }

    @Override
    public boolean equals(Object obj)
    {
        MapVertex vertex = (MapVertex)obj;
        double x1 = location.getX();
        double y1 = location.getY();
        double x2 = vertex.location.getX();
        double y2 = vertex.location.getY();
        boolean ret = (Double.compare(x1, x2) == 0) &&
                      (Double.compare(y1, y2) == 0);
        return ret;
    }

}

/**
 * This class represents a tuple of any vertex in the graph, and given
 * end vertex. We use this tuple as a key during graph search to
 * determine if we already have an optimal path from the current
 * vertex to the end vertex.
 * @author nishantsuneja
 *
 */
class PathLookupTuple {
    public PathLookupTuple(GeographicPoint curLocation,
                           GeographicPoint endLocation)
    {
        this.curLocation = curLocation;
        this.endLocation = endLocation;
    }

    public GeographicPoint getCurLocation()
    {
        return curLocation;
    }

    public GeographicPoint getEndLocation()
    {
        return endLocation;
    }

    @Override
    public int hashCode()
    {
        return new Double(23*curLocation.getX() + 29*curLocation.getY() +
               31*endLocation.getX() + 37*endLocation.getY()).intValue();
    }

    @Override
    public boolean equals(Object obj)
    {
        PathLookupTuple tuple = (PathLookupTuple)obj;
        double x1  = tuple.curLocation.getX();
        double y1 = tuple.curLocation.getY();
        double x2 = tuple.endLocation.getX();
        double y2 = tuple.endLocation.getY();

        double x3 = curLocation.getX();
        double y3 = curLocation.getY();
        double x4 = endLocation.getX();
        double y4 = endLocation.getY();

        boolean ret = (Double.compare(x1, x3) == 0) &&
                       (Double.compare(y1, y3) == 0) &&
                       (Double.compare(x2, x4) == 0) &&
                       (Double.compare(y2, y4) == 0);
        return ret;
    }

    @Override
    public String toString()
    {
        return "PathLookupTuple [curLocation=" + curLocation + ", endLocation="
                + endLocation + "]";
    }

    private final GeographicPoint curLocation;
    private final GeographicPoint endLocation;
}


/**
 * This class represents a graph of geographic locations.
 * Nodes in the graph are the intersections, and the edges are the road
 * segments.
 * @author nishantsuneja
 *
 */
public class MapGraph {

    private Map<GeographicPoint, MapVertex> vertexMap;
    // This set caches all the existing SHORTEST paths from a given
    // location to another given location in the map. We use this set
    // to optimizes the search algorithms, as we avoid running the
    // algorithm on the "common" part of the different searches, performed
    // on the same graph instance.
    private Set<PathLookupTuple> shortestPathExists;
    private int numVertices;
    private int numEdges;
    private enum SearchAlgo {DIJKSTRA, ASTAR};

    /**
     * This method resets the map graph to its earlier pristine state
     * before running a new search algorithm.
     */
    private void initializeGraph() {
        // Set the distances to inifinity.
        for (GeographicPoint location: vertexMap.keySet()) {
            MapVertex vertex = vertexMap.get(location);
            vertex.setDistanceFromOrigin(Double.MAX_VALUE);
        }
    }

    /**
     * This method validates the feasibility of adding an edge to the graph.
     * @return
     */
    private boolean validateEdge(GeographicPoint from,
                                 GeographicPoint to,
                                 String roadName,
                                 String roadType,
                                 double length)
    {
        if ((from == null) || (to == null)) { // Invalid endpoints
            return false;
        }

        // NOTE: The dataset seems to have empty road names too. So,
        // avoiding the check for empty roadnames.
        if (roadName == null) {
            return false; // Invalid road name.
        }

        // NOTE: The dataset seems to have empty road types too. So,
        // avoiding the check for empty roadtypes.
        if (roadType == null) {
            return false; // Invalid road type.
        }

        if (length <= 0) {  // Invalid edge length
            return false;
        }

        if (!vertexMap.containsKey(from) || !vertexMap.containsKey(to)) {
            return false; // One or both the endpoints are not present in the graph.
        }

        return true;
    }


    /**
     * This method returns the path between start and end vertices.
     * @param startVertex
     * @param endVertex
     * @param path
     */
    private void calculatePath(MapVertex startVertex,
                               MapVertex endVertex,
                               List<GeographicPoint> path)
    {
        MapVertex vertex;

        // Calculate the path.
        vertex = endVertex;
        while (vertex != null) {
            path.add(0, vertex.getLocation()); // Append to the beginning of the list.
            vertex = vertexMap.get(vertex.getParentLocation
                                   (endVertex.getLocation()));
        }
    }
 

    /** Find the path from start to goal using Dijkstra's algorithm/A* algorithm.
     * 
     * @param start The starting location
     * @param goal The goal location
     * @param queue The type of priority queue based upon the search
     *              algorithm.
     * @param nodeSearched A hook for visualization. 
     *        See assignment instructions for how to use it.
     * @return A tuple containing the following:
     *   a) The list of intersections that form the shortest path from 
     *      start to goal (including both start and goal).
     *   b) Number of nodes visited while determining the shortest path.
     */
    private Map.Entry<List<GeographicPoint>, Integer> findShortestPath(GeographicPoint start, 
                                                 GeographicPoint goal,
                                                PriorityQueue<MapVertex> queue,
                                                SearchAlgo searchAlgo,
                                                Consumer<GeographicPoint> nodeSearched)
    {
        Set<GeographicPoint> visitedSet = new HashSet<GeographicPoint>();
        List<GeographicPoint> path = new LinkedList<GeographicPoint>();
        MapVertex startVertex;
        MapVertex endVertex;
        MapVertex vertex; 
        double newDistanceFromOrigin;
        boolean pathFound;
        // Foe debugging only
        int numVerticesVisisted;

        // Lets first reset the map graph, before running the search,
        // because the user might run multiple searches on the same
        // graph instance.
        initializeGraph();

        // Ensure that the start and end points are present in the graph.
        if (!vertexMap.containsKey(start)) {
            return null;
        }
        startVertex = vertexMap.get(start);

        if (!vertexMap.containsKey(goal)) {
            return null;
        }
        endVertex = vertexMap.get(goal);

        // Begin with the source.
        vertex = startVertex;
        // Set the path distance to 0, as this is the starting vertex
        vertex.setDistanceFromOrigin(0.0);
        queue.add(vertex);

        /*
         * Expand the search after discovery of each new vertex, by choosing
         * the vertex which is closest to the source.
         */
        numVerticesVisisted = 0;
        pathFound = false;
        while (!queue.isEmpty()) {
            vertex = queue.remove();

            // For graph search visualization.
            nodeSearched.accept(vertex.getLocation());

            // For debugging only.
            numVerticesVisisted++;

            //System.out.println(searchAlgo.toString() + 
            //   " VISITING [NODE at location ("+ vertex.getLocation()+")]");

            if (!visitedSet.contains(vertex.getLocation())) {
                visitedSet.add(vertex.getLocation()); // Mark the vertex as visited.

                // We found the shortest path from start to this vertex.
                PathLookupTuple tuple1 = new PathLookupTuple(startVertex.getLocation(),
                                                            vertex.getLocation());
                shortestPathExists.add(tuple1);

                // Check if we already have a cached path from the current
                // vertex to the goal.
                PathLookupTuple tuple2 = new PathLookupTuple(vertex.getLocation(),
                                                    endVertex.getLocation());
                if (shortestPathExists.contains(tuple2)) {
                    pathFound = true;
                    System.out.println("Path from::" + vertex + " to " +
                            endVertex + " already cached. Terminating the search...");
                    break;
                }

                // Check if we reached our goal.
                if (vertex.getLocation().equals(endVertex.getLocation())) {
                    System.out.println("Found the vertex::" + endVertex +
                                       ". Terminating the search...");
                    pathFound = true;
                    break;
                }

                /*
                 * Iterate over all the edges of the most recently
                 * dequeued vertex.
                 */
                for (MapEdge edge: vertex.getOutEdges()) {
                    MapVertex neighbor = vertexMap.get(edge.getTo());
                    assert(neighbor != null);

                    // Let's not visit the neighbor which we have already
                    // visited, as this neighbor has already the shortest
                    // path set from the origin.
                    if (visitedSet.contains(neighbor.getLocation())) {
                        continue;
                    }

                    // Check if this neighbor can NOW be reached by a
                    // shorter route than before.
                    newDistanceFromOrigin = vertex.getDistanceFromOrigin() +
                            vertex.getNeighborEdgeLength(neighbor);
                    if (neighbor.getDistanceFromOrigin() > newDistanceFromOrigin) {
                        // Update the new distance from the source for
                        // the neighbor vertex.
                        neighbor.setDistanceFromOrigin(newDistanceFromOrigin);

                        // Update the parent for the neighbor, as we found
                        // a new route.
                        neighbor.setParentLocation(vertex.getLocation(),
                                                   endVertex.getLocation());

                        // Add the vertex's neighbor to the queue. 
                        queue.add(neighbor);
                    }
                }
            }
        }

        // If no path is found, return null.
        if (!pathFound) { // Path not found.
            System.out.println("No path found from::" + startVertex +
                    " to " + endVertex);
            return null;
        }

        // Calculate the path
        calculatePath(startVertex, endVertex, path);

        return new AbstractMap.SimpleEntry<List<GeographicPoint>, Integer>
                                               (path, numVerticesVisisted);
    }

    /**
     * For A star algorithm, the priority is the sum of distance
     * from source to a vertex and an "underestimation" of the
     * distance from the vertex to the destination. This ensures
     * faster convergence of the algorithm. To ensure the
     * underestimation, we use the straight line distance between
     * the 2 endpoints. 
     */
    private class AStarComparator implements Comparator<MapVertex> {

        public AStarComparator(GeographicPoint goal)
        {
            this.goal = goal;
        }

        @Override
        public int compare(MapVertex v1, MapVertex v2)
        {
             double dist1 = v1.getDistanceFromOrigin() +
                                v1.getLocation().distance(goal);
             double dist2 = v2.getDistanceFromOrigin() +
                                 v2.getLocation().distance(goal);

             if (dist1 < dist2) {
                 return -1;
             } else if (dist1 > dist2) {
                 return 1;
             } else {
                 return 0;
             }
        }

        private GeographicPoint goal; // The destination vertex.
    }

    /**
     * Comparator object to control the enqueue and dequeue operations into
     * our priority queue. The comparison between the 2
     * vertices is made based upon their distance from the origin. This
     * comparator method is used to calculate the shortest distance
     * when using dijakstra's algorithm.
     * @author nishantsuneja
     *
     */
    private class DjComparator implements Comparator<MapVertex> {

        @Override
        public int compare(MapVertex v1, MapVertex v2)
        {
             if (v1.getDistanceFromOrigin() < v2.getDistanceFromOrigin()) {
                 return -1;
             } else if (v1.getDistanceFromOrigin() > v2.getDistanceFromOrigin()) {
                 return 1;
             } else {
                 return 0;
             }
        }
    }


	/** 
	 * Create a new empty MapGraph 
	 */
	public MapGraph()
	{
		vertexMap = new HashMap<GeographicPoint, MapVertex>();
	    shortestPathExists = new LinkedHashSet<PathLookupTuple>();
		numVertices = 0;
		numEdges = 0;
	}
	
	/**
	 * Get the number of vertices (road intersections) in the graph
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices()
	{
	    return numVertices;
	}
	
	/**
	 * Return the intersections, which are the vertices in this graph.
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices()
	{
        return vertexMap.keySet();
	}
	
	/**
	 * Get the number of road segments in the graph
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges()
	{
        return numEdges;
	}

	/**Add a node corresponding to an intersection at a Geographic Point
	 * If the location is already in the graph or null, this method does 
	 * not change the graph.
	 * @param location  The location of the intersection
	 * @return true if a node was added, false if it was not (the node
	 * was already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location)
	{
	    if (location == null || vertexMap.containsKey(location)) {
	        return false;
	    }

	    // Initialize the vertex node.
	    MapVertex vertex = new MapVertex(location);

	    // Insert the new node into the vertex map.
	    vertexMap.put(location, vertex);

	    // Bump up the vertex count
	    numVertices++;

	    return true;
	}
	
	/**
	 * Adds a directed edge to the graph from pt1 to pt2.  
	 * Precondition: Both GeographicPoints have already been added to the graph
	 * @param from The starting point of the edge
	 * @param to The ending point of the edge
	 * @param roadName The name of the road
	 * @param roadType The type of the road
	 * @param length The length of the road, in km
	 * @throws IllegalArgumentException If the points have not already been
	 *   added as nodes to the graph, if any of the arguments is null,
	 *   or if the length is less than 0.
	 */
	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName,
			String roadType, double length) throws IllegalArgumentException {

	    if (!validateEdge(from, to, roadName, roadType, length)) {
	        throw new IllegalArgumentException("Invalid parameters while"
	                                   + " adding an edge to the graph.");
	    }

	    // First, lets create a directed edge between the 2 endpoints.
	    MapEdge edge = new MapEdge(from, to, roadName, roadType, length);

	    // Extract the start endpoint from the vertex map
	    MapVertex fromVertex = vertexMap.get(from);
	    assert(fromVertex != null);

	    // Hang the edge off the vertex
	    if (fromVertex.addOutEdge(edge)) {
	        // We successfully added an edge of the start endpoint. Bump
	        // up the edge count.
	        numEdges++;
	    }
	}


	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return bfs(start, goal, temp);
	}
	
	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, 
			 					     GeographicPoint goal,
			 					     Consumer<GeographicPoint> nodeSearched)
	{
	    Queue<MapVertex> bfsQueue = new LinkedList<MapVertex>();
	    Set<GeographicPoint> visitedSet = new HashSet<GeographicPoint>();
	    List<GeographicPoint> path = new LinkedList<GeographicPoint>();
	    MapVertex startVertex;
	    MapVertex endVertex;
	    MapVertex vertex;
	    boolean pathFound;

	    // Lets first reset the map graph, before running the search,
        // because the user might run multiple searches on the same
        // graph instance.
	    initializeGraph();

	    // Ensure that the start and end points are present in the graph.
	    if (!vertexMap.containsKey(start)) {
	        return null;
	    }
	    startVertex = vertexMap.get(start);

	    if (!vertexMap.containsKey(goal)) {
	        return null;
	    }
	    endVertex = vertexMap.get(goal);

	    // Begin with the source.
	    vertex = startVertex;
	    bfsQueue.add(vertex);

	    /*
	     * Move breath wise through the graph till we encounter the
	     * goal vertex.
	     */
	    pathFound = false;
	    while (!bfsQueue.isEmpty()) {
	        vertex = bfsQueue.remove();
	        // For visualization of the graph traversal
	        nodeSearched.accept(vertex.getLocation());
            visitedSet.add(vertex.getLocation()); // Mark the vertex as visited.
 
            // We found the shortest path from start to this vertex.
            shortestPathExists.add(new PathLookupTuple(startVertex.getLocation(),
                                               vertex.getLocation()));

	        // Check if we reached our goal.
	        if (vertex.getLocation().equals(endVertex.getLocation()) ||
                    // We already have a cached path from the current
                    // vertex to the goal.
                    shortestPathExists.contains(new PathLookupTuple
                                        (vertex.getLocation(),
                                        endVertex.getLocation()))) {
	            pathFound = true;
	            break;
	        }

	        /*
	         * Iterate over all the edges of the most recently
	         * dequeued vertex.
	         */
	        for (MapEdge edge: vertex.getOutEdges()) {
	            // Iterate over the neighboring vertex only if we didn't
	            // already visit it.
	            if (!visitedSet.contains(edge.getTo())) {
	                MapVertex neighbor = vertexMap.get(edge.getTo());
	                assert(neighbor != null);

	                // Enqueue the newly discovered vertex.
	                bfsQueue.add(neighbor);
	                // Maintain the parent child relationship
	                neighbor.setParentLocation(vertex.getLocation(),
	                                           endVertex.getLocation());
	            }
	        }
	    }

	    // If no path is found, return null.
	    if (!pathFound) { // Path not found.
	        System.out.println("No path found from::" + startVertex +
	                           " to " + endVertex);
	        return null;
	    }

	    // Calculate the path from start to end.
	    calculatePath(startVertex, endVertex, path);

		return path;
	}


	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		// You do not need to change this method.
        Consumer<GeographicPoint> temp = (x) -> {};
        return dijkstra(start, goal, temp);
	}


    /** Find the path from start to goal using Dijkstra's algorithm
     * 
     * @param start The starting location
     * @param goal The goal location
     * @return The list of intersections that form the shortest path from 
     *   start to goal (including both start and goal).
     */	
    public List<GeographicPoint> dijkstra(GeographicPoint start, 
                                          GeographicPoint goal,
                                   Consumer<GeographicPoint> nodeSearched)
    {
         DjComparator djComp = new DjComparator();
         PriorityQueue<MapVertex> djQueue = new PriorityQueue<MapVertex>
                  (getNumVertices(), djComp);
         Map.Entry<List<GeographicPoint>, Integer> path =
                 findShortestPath(start, goal, djQueue,
                              SearchAlgo.DIJKSTRA, nodeSearched);
         System.out.println("Number of nodes visisted via DIJKSTRA algo::"
                                                       + path.getValue());
         return path.getKey();
    }
	
	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return aStarSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, 
											 GeographicPoint goal,
											 Consumer<GeographicPoint> nodeSearched)
	{
	    AStarComparator aStarComp = new AStarComparator(goal);
        PriorityQueue<MapVertex> aStarQueue = new PriorityQueue<MapVertex>
                                              (getNumVertices(), aStarComp);
        Map.Entry<List<GeographicPoint>, Integer> path =
                findShortestPath(start, goal, aStarQueue,
                                SearchAlgo.ASTAR, nodeSearched);
        System.out.println("Number of nodes visited via A* algo::"
                                                       + path.getValue());
        return path.getKey();
	}


	public static void main(String[] args)
	{
	    /*Set<MapVertex> gSet = new HashSet<MapVertex>();
	    Set<PathLookupTuple> tSet = new HashSet<PathLookupTuple>();
        GeographicPoint g1 = new GeographicPoint(8.0, -1.0);
        GeographicPoint g2 = new GeographicPoint(8.0, -1.0);

        MapVertex v1 = new MapVertex(g1);
        MapVertex v2 = new MapVertex(g2);

        System.out.println(v1.hashCode() + "," + v2.hashCode());
        System.out.println(v1.equals(v2));

        gSet.add(v1);
        gSet.add(v2);
        System.out.println(gSet);

        PathLookupTuple t1 = new PathLookupTuple(g1, g2);
        PathLookupTuple t2 = new PathLookupTuple(g1, g2);

        tSet.add(t1);
        tSet.add(t2);
        System.out.println(tSet);*/

	    MapGraph simpleTestMap = new MapGraph();
        GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);

        GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
        GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);

        System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
        List<GeographicPoint> testroute = simpleTestMap.dijkstra(testStart,testEnd);
        List<GeographicPoint> testroute2 = simpleTestMap.aStarSearch(testStart,testEnd);
        System.out.println(testroute);
        System.out.println(testroute2);
        assert(testroute.equals(testroute2));

        MapGraph testMap = new MapGraph();
        GraphLoader.loadRoadMap("data/maps/utc.map", testMap);

        // A very simple test using real data
        testStart = new GeographicPoint(32.869423, -117.220917);
        testEnd = new GeographicPoint(32.869255, -117.216927);
        System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
        testroute = testMap.dijkstra(testStart,testEnd);
        testroute2 = testMap.aStarSearch(testStart,testEnd);
        System.out.println(testroute);
        System.out.println(testroute2);

        // A slightly more complex test using real data
        testStart = new GeographicPoint(32.8674388, -117.2190213);
        testEnd = new GeographicPoint(32.8697828, -117.2244506);
        System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
        testroute = testMap.dijkstra(testStart,testEnd);
        testroute2 = testMap.aStarSearch(testStart,testEnd);
        System.out.println(testroute);
        System.out.println(testroute2);

	}
	
}

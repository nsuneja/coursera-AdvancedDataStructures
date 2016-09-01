/**
 * @author UCSD MOOC development team and YOU
 * This class represents a graph of geographic locations.
 * Nodes in the graph are the intersections, and the edges are the road
 * segments.
 * 
 */
package roadgraph;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
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
    private GeographicPoint location;
    private GeographicPoint parentLocation;
    private List<MapEdge> outEdges;

    public MapVertex(GeographicPoint location) {
        this.location = location;
        this.parentLocation = null;
        this.outEdges = new ArrayList<MapEdge>();
    }

    public GeographicPoint getLocation()
    {
        return location;
    }
    public void setLocation(GeographicPoint location)
    {
        this.location = location;
    }
    public List<MapEdge> getOutEdges()
    {
        return outEdges;
    }
    public void setOutEdges(List<MapEdge> outEdges)
    {
        this.outEdges = outEdges;
    }

    public GeographicPoint getParentLocation()
    {
        return parentLocation;
    }

    public void setParentLocation(GeographicPoint parentLocation)
    {
        this.parentLocation = parentLocation;
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
    public boolean equals(Object obj)
    {
        MapVertex vertex = (MapVertex)obj;
        return location.equals(vertex.getLocation());
    }

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
    private int numVertices;
    private int numEdges;

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
	 * Create a new empty MapGraph 
	 */
	public MapGraph()
	{
		vertexMap = new HashMap<GeographicPoint, MapVertex>();
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
			 					     GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
	    Queue<MapVertex> bfsQueue = new LinkedList<MapVertex>();
	    Set<GeographicPoint> visitedSet = new HashSet<GeographicPoint>();
	    List<GeographicPoint> path = new LinkedList<GeographicPoint>();
	    MapVertex startVertex;
	    MapVertex endVertex;
	    MapVertex vertex; 

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
	    while (!bfsQueue.isEmpty()) {
	        vertex = bfsQueue.remove();
	        nodeSearched.accept(vertex.getLocation());
            visitedSet.add(vertex.getLocation()); // Mark the vertex as visited.

	        // Check if we reached our goal.
	        if (vertex.getLocation().equals(goal)) {
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
	                neighbor.setParentLocation(vertex.getLocation());
	            }
	        }
	    }

	    // If no path is found, return null.
	    if (!vertex.equals(endVertex)) { // Path not found.
	        return null;
	    }

	    // Calculate the path.
	    assert(vertex.equals(endVertex)); // Safety check.
	    while (vertex.getParentLocation() != null) {
	        path.add(0, vertex.getLocation()); // Append to the beginning of the list.
	        vertex = vertexMap.get(vertex.getParentLocation());
	    }
	
	    // We reached the start point.
	    assert(vertex.equals(startVertex));
	    path.add(0, startVertex.getLocation());
	
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
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, 
										  GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// TODO: Implement this method in WEEK 3

		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());
		
		return null;
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
											 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// TODO: Implement this method in WEEK 3
		
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());
		
		return null;
	}

	
	
	public static void main(String[] args)
	{
		System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", firstMap);
		System.out.println("DONE.");
		
		// You can use this method for testing.  
		
		
		/* Here are some test cases you should try before you attempt 
		 * the Week 3 End of Week Quiz, EVEN IF you score 100% on the 
		 * programming assignment.
		 */
		/*
		MapGraph simpleTestMap = new MapGraph();
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);
		
		GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
		GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);
		
		System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
		List<GeographicPoint> testroute = simpleTestMap.dijkstra(testStart,testEnd);
		List<GeographicPoint> testroute2 = simpleTestMap.aStarSearch(testStart,testEnd);
		
		
		MapGraph testMap = new MapGraph();
		GraphLoader.loadRoadMap("data/maps/utc.map", testMap);
		
		// A very simple test using real data
		testStart = new GeographicPoint(32.869423, -117.220917);
		testEnd = new GeographicPoint(32.869255, -117.216927);
		System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
		testroute = testMap.dijkstra(testStart,testEnd);
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		
		
		// A slightly more complex test using real data
		testStart = new GeographicPoint(32.8674388, -117.2190213);
		testEnd = new GeographicPoint(32.8697828, -117.2244506);
		System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
		testroute = testMap.dijkstra(testStart,testEnd);
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		*/
		
		
		/* Use this code in Week 3 End of Week Quiz */
		/*MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");

		GeographicPoint start = new GeographicPoint(32.8648772, -117.2254046);
		GeographicPoint end = new GeographicPoint(32.8660691, -117.217393);
		
		
		List<GeographicPoint> route = theMap.dijkstra(start,end);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);

		*/
		
	}
	
}

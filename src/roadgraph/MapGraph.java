/**
 * @author UCSD MOOC development team and YOU
 * 
 * A class which reprsents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
package roadgraph;


import java.util.*;
import java.util.function.Consumer;

import geography.GeographicPoint;
import util.GraphLoader;

/**
 * @author UCSD MOOC development team and eliyooy
 * 
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
public class MapGraph {
	//TODO: Add your member variables here in WEEK 2

	private int numVertices;
	private int numEdges;
	private Set<GeographicPoint> vertices;
	private HashMap<GeographicPoint, ArrayList<MapEdge>> edges;

	/** 
	 * Create a new empty MapGraph
	 * numVertices - The number of vertices on the map
	 * numEdges - The number of edges in the map
	 * vertices - Set containing all vertices in the map
	 * edges - HashMap containing all edges in the map, keyed to the origin GeographicPoint of the edge.
	 * Key returns an ArrayList with all edges extending from the keyed GeographicPoint.
	 * Each edge is represented by a MapEdge object which contains all data for the edge.
	 */
	public MapGraph()
	{
		// TODO: Implement in this constructor in WEEK 2
		numVertices = 0;
		numEdges = 0;
		vertices = new HashSet<>();
		edges = new HashMap<>();
	}
	
	/**
	 * Get the number of vertices (road intersections) in the graph
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices()
	{
		//TODO: Implement this method in WEEK 2
		return numVertices;
	}
	
	/**
	 * Return the intersections, which are the vertices in this graph.
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices()
	{
		//TODO: Implement this method in WEEK 2
		return vertices;
	}
	
	/**
	 * Get the number of road segments in the graph
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges()
	{
		//TODO: Implement this method in WEEK 2
		return numEdges;
	}

	
	
	/** Add a node corresponding to an intersection at a Geographic Point
	 * If the location is already in the graph or null, this method does 
	 * not change the graph.
	 * @param location  The location of the intersection
	 * @return true if a node was added, false if it was not (the node
	 * was already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location) {
		// TODO: Implement this method in WEEK 2
		if (graphContainsNode(location) || location == null) {
			return false;
		} else {
			vertices.add(location);
			numVertices += 1;
			return true;
		}
	}

	public boolean graphContainsNode(GeographicPoint location) {

		if(vertices.contains(location)) {
			return true;
		}

		return false;
	}
	
	/**
	 * Adds a directed edge to the graph from pt1 to pt2.  
	 * Precondition: Both GeographicPoints have already been added to the graph
	 * Creates a new MapEdge object with param data. Adds the new edge to the HashMap
	 * of edges.
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

		//TODO: Implement this method in WEEK 2
		if(!graphContainsNode(from) || !graphContainsNode(to) || length < 0 || from == null || to == null ||
				roadName == null || roadType == null) {
			throw new IllegalArgumentException("Please enter valid values for all parameters.");
		}

		MapEdge newEdge = new MapEdge(from, to, roadName, roadType, length);

		if(edges.containsKey(from)) {
			edges.get(from).add(newEdge);
		} else {
			edges.put(from, new ArrayList<>());
			edges.get(from).add(newEdge);
		}

		numEdges += 1;
		
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
		// TODO: Implement this method in WEEK 2
		
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next);

		Queue<GeographicPoint> q = new LinkedList<>();
		LinkedList<GeographicPoint> visited = new LinkedList<>();
		List<ArrayList<GeographicPoint>> paths = new LinkedList<>();
		loadPaths(paths, q, visited, start);
		int possiblePaths = 1;

		while(!q.isEmpty()) {
			if (!edges.containsKey(q.peek())) {
				q.remove();
				continue;
			}

			GeographicPoint curr = q.remove();
			nodeSearched.accept(curr);

			if (curr.equals(goal)) {

				return determineShortestList(paths, curr, goal);

			} else {

				for (MapEdge neighbor : edges.get(curr)) {
					if (!visited.contains(neighbor.getTo())) {
						visited.add(neighbor.getTo());
						q.add(neighbor.getTo());
					}

					possiblePaths += processNewPossiblePaths(paths, curr, neighbor, possiblePaths);

				}
			}
		}
		return null;
	}

	public List<ArrayList<GeographicPoint>> loadPaths(List<ArrayList<GeographicPoint>> paths,
													 Queue<GeographicPoint> q, LinkedList<GeographicPoint> visited,
													  GeographicPoint start) {
		ArrayList<GeographicPoint> firstList = new ArrayList<>();
		firstList.add(start);
		paths.add(firstList);
		q.add(start);
		visited.add(start);

		return paths;
	}

	public List<ArrayList<GeographicPoint>> determinePossiblePaths(List<ArrayList<GeographicPoint>> paths,
																   GeographicPoint curr) {
		List<ArrayList<GeographicPoint>> possiblePaths = new LinkedList<>();

		for (ArrayList<GeographicPoint> currList : paths) {
			GeographicPoint lastPoint = currList.get(currList.size() - 1);

			if(edges.get(lastPoint) == null) {
				continue;
			}

			for (MapEdge currEdge : edges.get(lastPoint)) {
				if (currEdge.getTo().equals(curr)) {
					possiblePaths.add(currList);
					break;
				}
			}
		}

		return possiblePaths;
	}

	public ArrayList<GeographicPoint> determineShortestList(List<ArrayList<GeographicPoint>> paths,
															GeographicPoint curr,
															GeographicPoint goal) {
		List<ArrayList<GeographicPoint>> possiblePaths = determinePossiblePaths(paths, curr);
		ArrayList<GeographicPoint> shortestList = possiblePaths.get(0);

		for (ArrayList<GeographicPoint> finalist : possiblePaths) {
			if (finalist.size() < shortestList.size()) {
				shortestList = finalist;
			}
		}

		shortestList.add(goal);
		return shortestList;
	}

	public int processNewPossiblePaths(List<ArrayList<GeographicPoint>> paths, GeographicPoint curr, MapEdge neighbor,
									 int possiblePaths) {
		int possibleNewPaths = 0;

		for (int i=0; i<possiblePaths; i++) {
			if (paths.get(i).get(paths.get(i).size() - 1).equals(curr)) {
				ArrayList<GeographicPoint> newList = new ArrayList<>();
				newList.addAll(paths.get(i));
				newList.add(neighbor.getTo());
				paths.add(newList);
				possibleNewPaths += 1;
			}
		}

		return possibleNewPaths;

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
		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", theMap);
		System.out.println("DONE.");
		
		// You can use this method for testing.
		GeographicPoint start = new GeographicPoint(4.0, 2.0);
		GeographicPoint end = new GeographicPoint(8.0, -1.0);
		theMap.bfs(start, end);

		/* Use this code in Week 3 End of Week Quiz
		MapGraph theMap = new MapGraph();
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

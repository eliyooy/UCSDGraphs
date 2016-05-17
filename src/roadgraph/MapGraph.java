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
	private ArrayList<GeographicPoint> recordVisitedNodes = new ArrayList<>();

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
		loadPathsBFS(paths, q, visited, start);
		int possiblePaths = 1;

		while(!q.isEmpty()) {
			if (!edges.containsKey(q.peek())) {
				q.remove();
				continue;
			}

			GeographicPoint curr = q.remove();
			nodeSearched.accept(curr);
			recordVisitedNodes.add(curr);

			if (curr.equals(goal)) {

				return determineShortestListBFS(paths, curr, goal);

			} else {

				for (MapEdge neighbor : edges.get(curr)) {
					if (!visited.contains(neighbor.getTo())) {
						visited.add(neighbor.getTo());
						q.add(neighbor.getTo());
					}

					possiblePaths += processNewPossiblePathsBFS(paths, curr, neighbor, possiblePaths);

				}
			}
		}
		return null;
	}

	public void loadPathsBFS(List<ArrayList<GeographicPoint>> paths, Queue<GeographicPoint> q,
						  LinkedList<GeographicPoint> visited, GeographicPoint start) {
		ArrayList<GeographicPoint> firstList = new ArrayList<>();
		firstList.add(start);
		paths.add(firstList);
		q.add(start);
		visited.add(start);
	}

	public int processNewPossiblePathsBFS(List<ArrayList<GeographicPoint>> paths, GeographicPoint curr, MapEdge neighbor,
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

	public ArrayList<GeographicPoint> determineShortestListBFS(List<ArrayList<GeographicPoint>> paths,
															GeographicPoint curr,
															GeographicPoint goal) {
		List<ArrayList<GeographicPoint>> possiblePaths = determinePossiblePathsAllSearches(paths, curr);
		ArrayList<GeographicPoint> shortestList = possiblePaths.get(0);

		for (ArrayList<GeographicPoint> finalist : possiblePaths) {
			if (finalist.size() < shortestList.size()) {
				shortestList = finalist;
			}
		}

		shortestList.add(goal);
		return shortestList;
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

		LinkedList<GeographicPoint> q = new LinkedList<>();
		LinkedList<GeographicPoint> visited = new LinkedList<>();
		HashMap<ArrayList<GeographicPoint>, Double> paths = new HashMap<>();
		List<ArrayList<GeographicPoint>> pathKeys = new LinkedList<>();
		loadPathsDijkstraAStar(paths, pathKeys, q, visited, start);
		int possiblePaths = 1;

		while(!q.isEmpty()) {
			if (!edges.containsKey(q.peek())) {
				q.remove();
				continue;
			}

			GeographicPoint curr = q.removeFirst();
			nodeSearched.accept(curr);
			recordVisitedNodes.add(curr);

			if (curr.getX() == goal.getX() && curr.getY() == goal.getY()) {

				return determineShortestListDijkstraAStar(paths, pathKeys, curr, goal);

			} else {

				for (MapEdge neighbor : edges.get(curr)) {
					if (!visited.contains(neighbor.getTo())) {
						visited.add(neighbor.getTo());
						q.add(neighbor.getTo());
						Collections.sort(q, (a, b) -> a.distance(start) < b.distance(start) ? -1 :
								a.distance(start) == b.distance(start) ? 0 : 1);
					}

					possiblePaths += processNewPossiblePathsDijkstraAStar(paths, pathKeys, curr, neighbor, possiblePaths);

				}
			}
		}
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
		LinkedList<GeographicPoint> q = new LinkedList<>();
		LinkedList<GeographicPoint> visited = new LinkedList<>();
		HashMap<ArrayList<GeographicPoint>, Double> paths = new HashMap<>();
		List<ArrayList<GeographicPoint>> pathKeys = new LinkedList<>();
		loadPathsDijkstraAStar(paths, pathKeys, q, visited, start);
		int possiblePaths = 1;

		while(!q.isEmpty()) {
			if (!edges.containsKey(q.peek())) {
				q.remove();
			}

			GeographicPoint curr = q.removeFirst();
			nodeSearched.accept(curr);
			recordVisitedNodes.add(curr);

			if (curr.getX() == goal.getX() && curr.getY() == goal.getY()) {
				return determineShortestListDijkstraAStar(paths, pathKeys, curr, goal);

			} else {

				for (MapEdge neighbor : edges.get(curr)) {
					if (!visited.contains(neighbor.getTo())) {
						visited.add(neighbor.getTo());
						q.add(neighbor.getTo());
						Collections.sort(q, (a, b) -> (a.distance(start) + a.distance(goal)) <
								(b.distance(start) + b.distance(goal)) ? -1 :
								(a.distance(start) + a.distance(goal)) == (b.distance(start) + b.distance(goal)) ? 0 : 1);
					}

					possiblePaths += processNewPossiblePathsDijkstraAStar(paths, pathKeys, curr, neighbor, possiblePaths);

				}
			}
		}
		return null;
	}

	public void loadPathsDijkstraAStar(HashMap<ArrayList<GeographicPoint>, Double> paths,
									   List<ArrayList<GeographicPoint>> pathKeys, Queue<GeographicPoint> q,
									   LinkedList<GeographicPoint> visited, GeographicPoint start) {
		ArrayList<GeographicPoint> firstList = new ArrayList<>();
		firstList.add(start);
		pathKeys.add(firstList);
		paths.put(firstList, 0.0);
		q.add(start);
		visited.add(start);
	}


	public int processNewPossiblePathsDijkstraAStar(HashMap<ArrayList<GeographicPoint>, Double> paths,
													List<ArrayList<GeographicPoint>> pathKeys,
													GeographicPoint curr, MapEdge neighbor, int possiblePaths) {
		int possibleNewPaths = 0;

		for (int i=0; i<possiblePaths; i++) {
			if (pathKeys.get(i).get(pathKeys.get(i).size() - 1).equals(curr)) {
				ArrayList<GeographicPoint> newList = new ArrayList<>();
				newList.addAll(pathKeys.get(i));
				newList.add(neighbor.getTo());
				pathKeys.add(newList);
				paths.put(newList, paths.get(pathKeys.get(i)) + neighbor.getLength());
				possibleNewPaths += 1;
			}
		}

		return possibleNewPaths;
	}

	public ArrayList<GeographicPoint> determineShortestListDijkstraAStar(HashMap<ArrayList<GeographicPoint>, Double> paths,
																		 List<ArrayList<GeographicPoint>> pathKeys,
																		 GeographicPoint curr,
																		 GeographicPoint goal) {
		List<ArrayList<GeographicPoint>> possiblePaths = determinePossiblePathsAllSearches(pathKeys, curr);
		ArrayList<GeographicPoint> shortestList = possiblePaths.get(0);

		for(ArrayList<GeographicPoint> currList : possiblePaths) {
			Double newDistance = paths.get(currList) + currList.get(currList.size() - 1).distance(goal);
			currList.add(goal);

			for(int i=1; i<currList.size(); i++) {
				newDistance += currList.get(i).distance(currList.get(i-1));
			}

			paths.put(currList, newDistance);
		}

		for (ArrayList<GeographicPoint> finalist : possiblePaths) {
			if (paths.get(finalist) < paths.get(shortestList)) {
				shortestList = finalist;
			}
		}

		return shortestList;
	}

	public List<ArrayList<GeographicPoint>> determinePossiblePathsAllSearches(List<ArrayList<GeographicPoint>> paths,
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

	public void tightSearchScoreGenerator(GeographicPoint start, GeographicPoint goal) {

		double BFSTightSearchScore = processTightSearchScore(start, goal, (ArrayList<GeographicPoint>) bfs(start, goal));
		double dijkstraTightSearchScore = processTightSearchScore(start, goal, (ArrayList<GeographicPoint>) dijkstra(start, goal));
		double aStarTightSearchScore = processTightSearchScore(start, goal, (ArrayList<GeographicPoint>) aStarSearch(start, goal));

		System.out.println("Average Tight Search Scores: \nBFS: " + BFSTightSearchScore + "\nDijkstra: " +
		dijkstraTightSearchScore + "\nAStar: " + aStarTightSearchScore);

	}

	public double processTightSearchScore(GeographicPoint start, GeographicPoint goal, ArrayList<GeographicPoint> route) {
		GeographicPoint closestPoint = start;
		double totalPointDistances = 0.0;
		for( GeographicPoint currPoint : recordVisitedNodes ) {
			for( GeographicPoint routePoint : route ) {
				if( currPoint.distance(routePoint) < currPoint.distance(closestPoint) ) {
					closestPoint = routePoint;
				}
			}

			totalPointDistances += currPoint.distance(closestPoint);
		}

		double tightSearchScore = (totalPointDistances / recordVisitedNodes.size());

		recordVisitedNodes.clear();

		return tightSearchScore;
	}

	public static void main(String[] args)
	{
		System.out.print("Making a new map...");
		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/san_diego.map", theMap);
		System.out.println("DONE.");
		
		/*// You can use this method for testing.
		GeographicPoint start = new GeographicPoint(4.0, 2.0);
		GeographicPoint end = new GeographicPoint(8.0, -1.0);
		theMap.bfs(start, end);

		// Use this code in Week 3 End of Week Quiz
		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");*/

		GeographicPoint start = new GeographicPoint(32.7229776, -117.1665031);
		GeographicPoint end = new GeographicPoint(32.717812, -117.1601757);
		
		
		List<GeographicPoint> route = theMap.dijkstra(start,end);
		System.out.println(route);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);
		System.out.println(route2);

		theMap.tightSearchScoreGenerator(start, end);


		
	}
	
}

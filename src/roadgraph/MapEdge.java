package roadgraph;

import geography.GeographicPoint;

/**
 * Created by eliyahudorsky on 5/3/16.
 */
public class MapEdge {

    private GeographicPoint from;
    private GeographicPoint to;
    private String roadName;
    private String roadType;
    private double length;

    public MapEdge(GeographicPoint from, GeographicPoint to, String roadName, String roadType, double length) {
        this.from = from;
        this.to = to;
        this.roadName = roadName;
        this.roadType = roadType;
        this.length = from.distance(to);
    }

    public GeographicPoint getFrom() {
        return from;
    }

    public void setFrom(GeographicPoint from) {
        this.from = from;
    }

    public GeographicPoint getTo() {
        return to;
    }

    public void setTo(GeographicPoint to) {
        this.to = to;
    }

    public String getRoadName() {
        return roadName;
    }

    public void setRoadName(String roadName) {
        this.roadName = roadName;
    }

    public String getRoadType() {
        return roadType;
    }

    public void setRoadType(String roadType) {
        this.roadType = roadType;
    }

    public double getLength() {
        return length;
    }

    public void setLength(double length) {
        this.length = length;
    }
}

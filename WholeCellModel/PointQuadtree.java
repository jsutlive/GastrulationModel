import java.util.ArrayList;
import java.util.List;

/**
 * A point quadtree: stores an element at a 2D position, 
 * with children at the subdivided quadrants
 * 
 * @author Shicheng Huang, F002BMN, CS10
 * 
 * @author Chris Bailey-Kellogg, Dartmouth CS 10, Spring 2015
 * @author CBK, Spring 2016, explicit rectangle
 * @author CBK, Fall 2016, generic with Point2D interface
 * 
 */
public class PointQuadtree<E extends Point2D> {
	private E point;							// the point anchoring this node
	private int x1, y1;							// upper-left corner of the region
	private int x2, y2;							// bottom-right corner of the region
	private PointQuadtree<E> c1, c2, c3, c4;	// children

	/**
	 * Initializes a leaf quadtree, holding the point in the rectangle
	 */
	public PointQuadtree(E point, int x1, int y1, int x2, int y2) {
		this.point = point;
		this.x1 = x1; this.y1 = y1; this.x2 = x2; this.y2 = y2;
	}

	// Getters
	
	public E getPoint() {
		return point;
	}

	public int getX1() {
		return x1;
	}

	public int getY1() {
		return y1;
	}

	public int getX2() {
		return x2;
	}

	public int getY2() {
		return y2;
	}

	/**
	 * Returns the child (if any) at the given quadrant, 1-4
	 * @param quadrant	1 through 4
	 */
	public PointQuadtree<E> getChild(int quadrant) {
		if (quadrant==1) return c1;
		if (quadrant==2) return c2;
		if (quadrant==3) return c3;
		if (quadrant==4) return c4;
		return null;
	}

	/**
	 * Returns whether or not there is a child at the given quadrant, 1-4
	 * @param quadrant	1 through 4
	 */
	public boolean hasChild(int quadrant) {
		return (quadrant==1 && c1!=null) || (quadrant==2 && c2!=null) || (quadrant==3 && c3!=null) || (quadrant==4 && c4!=null);
	}

	/**
	 * Inserts the point into the tree
	 */
	public void insert(E p2) {
		// TODO: YOUR CODE HERE
		if(p2.getX()>point.getX()) {
			if(p2.getY()>point.getY()) {
				if(getChild(4)==null) {
					c4 = new PointQuadtree<E>(p2,(int)(point.getX()),(int)(point.getY()),x2,y2);
				}else {
					c4.insert(p2);
				}
			}else {
				if(getChild(1)==null) {
					c1 = new PointQuadtree<E>(p2,(int)(point.getX()),y1,x2,(int)(point.getY()));
				}else {
					c1.insert(p2);
				}
			}
		}else {
			if(p2.getY()>point.getY()) {
				if(getChild(3)==null) {
					c3 = new PointQuadtree<E>(p2,x1,(int)(point.getY()),(int)(point.getX()),y2);
				}else {
					c3.insert(p2);
				}	
			}else {
				if(getChild(2)==null) {
					c2 = new PointQuadtree<E>(p2,x1,y1,(int)(point.getX()),(int)(point.getY()));
				}else {
					c2.insert(p2);
				}
			}
		}
	}
	
	/**
	 * Finds the number of points in the quadtree (including its descendants)
	 */
	public int size() {
		// TODO: YOUR CODE HERE
		int num = 1;
		if(getChild(1)!=null) num+=c1.size();
		if(getChild(2)!=null) num+=c2.size();
		if(getChild(3)!=null) num+=c3.size();
		if(getChild(4)!=null) num+=c4.size();
		return num;
	}
	
	/**
	 * Builds a list of all the points in the quadtree (including its descendants)
	 */
	public List<E> allPoints() {
		// TODO: YOUR CODE HERE
		List<E> points = new ArrayList<E>();
		allPointsHelper(points);
		return points;
	}	
	
	public void allPointsHelper(List<E> points) {
		points.add(point);
		if(getChild(1)!=null) c1.allPointsHelper(points);
		if(getChild(2)!=null) c2.allPointsHelper(points);
		if(getChild(3)!=null) c3.allPointsHelper(points);
		if(getChild(4)!=null) c4.allPointsHelper(points);
		
	}

	/**
	 * Uses the quadtree to find all points within the circle
	 * @param cx	circle center x
	 * @param cy  	circle center y
	 * @param cr  	circle radius
	 * @return    	the points in the circle (and the qt's rectangle)
	 */
	public List<E> findInCircle(double cx, double cy, double cr) {
		// TODO: YOUR CODE HERE
		List<E> points = new ArrayList<E>();
		findInCircleHelper(points,cx,cy,cr);
		return points;
	}
	
	public void findInCircleHelper(List<E> points, double cx, double cy, double cr) {
		if(Geometry.circleIntersectsRectangle(cx, cy, cr, (double)x1, (double)y1, (double)x2, (double)y2)) {
			if(Geometry.pointInCircle(point.getX(), point.getY(), cx, cy, cr)) points.add(point);
			if(getChild(1)!=null) c1.findInCircleHelper(points,cx,cy,cr);
			if(getChild(2)!=null) c2.findInCircleHelper(points,cx,cy,cr);
			if(getChild(3)!=null) c3.findInCircleHelper(points,cx,cy,cr);
			if(getChild(4)!=null) c4.findInCircleHelper(points,cx,cy,cr);
		}
	}

	// TODO: YOUR CODE HERE for any helper methods
}

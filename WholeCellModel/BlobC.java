import java.awt.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;


public class BlobC extends Blob implements Point2D{
	public DecimalFormat df;
	private double ax = 0, ay = 0;
	private static double arrowScale = 300;
	private ArrayList<Double> xLoad;
	private ArrayList<Double> yLoad;
	private ArrayList<Double> xLoadBlue;
	private ArrayList<Double> yLoadBlue;
	private ArrayList<Double> xLoadGreen;
	private ArrayList<Double> yLoadGreen;
	private ArrayList<Double> xLoadCyan;
	private ArrayList<Double> yLoadCyan;
	
	public static void drawArrows(boolean drwAws,double arrScl) {
		arrowScale = arrScl;
	}
	
	public BlobC() {
		super();
		df = new DecimalFormat("#.##");
		df.setRoundingMode(RoundingMode.DOWN);
		xLoad = new ArrayList<Double>();
		yLoad = new ArrayList<Double>();
		xLoadBlue = new ArrayList<Double>();
		yLoadBlue = new ArrayList<Double>();
		xLoadGreen = new ArrayList<Double>();
		yLoadGreen = new ArrayList<Double>();
		xLoadCyan = new ArrayList<Double>();
		yLoadCyan = new ArrayList<Double>();
	}
	public BlobC(double x, double y) {
		super(x,y);
		df = new DecimalFormat("#.##");
		df.setRoundingMode(RoundingMode.DOWN);
		xLoad = new ArrayList<Double>();
		yLoad = new ArrayList<Double>();
		xLoadBlue = new ArrayList<Double>();
		yLoadBlue = new ArrayList<Double>();
		xLoadGreen = new ArrayList<Double>();
		yLoadGreen = new ArrayList<Double>();
		xLoadCyan = new ArrayList<Double>();
		yLoadCyan = new ArrayList<Double>();
	}
	public BlobC(double x, double y, double r) {
		super(x,y,r);
		df = new DecimalFormat("#.##");
		df.setRoundingMode(RoundingMode.DOWN);
		xLoad = new ArrayList<Double>();
		yLoad = new ArrayList<Double>();
		xLoadBlue = new ArrayList<Double>();
		yLoadBlue = new ArrayList<Double>();
		xLoadGreen = new ArrayList<Double>();
		yLoadGreen = new ArrayList<Double>();
		xLoadCyan = new ArrayList<Double>();
		yLoadCyan = new ArrayList<Double>();
	}
	
	public double getVelocityX() {
		return dx;
	}
	
	public double getVelocityY() {
		return dy;
	}
	
	public double getAccelerationX() {
		return ax;
	}
	
	public double getAccelerationY() {
		return ay;
	}
	
	public void setAcceleration(double ax, double ay) {
		this.ax = ax;
		this.ay = ay;
	}
	
	public void applyLoad(double l1, double l2, boolean drawA) {
//		l1 = Math.floor(l1*100)/100;
//		l2 = Math.floor(l2*100)/100;
		ax += l1;
		ay += l2;
		if(drawA) {
			xLoad.add(l1);
			yLoad.add(l2);
		}
	}
	
	public void applyLoad(double l1, double l2, Color color) {
//		l1 = Math.floor(l1*100)/100;
//		l2 = Math.floor(l2*100)/100;
		ax += l1;
		ay += l2;
		if(color == Color.BLUE) {
			xLoadBlue.add(l1);
			yLoadBlue.add(l2);
		}
		else if(color == Color.GREEN) {
			xLoadGreen.add(l1);
			yLoadGreen.add(l2);
		}
		else if(color == Color.CYAN) {
			xLoadCyan.add(l1);
			yLoadCyan.add(l2);
		}
		else {
			xLoad.add(l1);
			yLoad.add(l2);
		}
	}
	
	public void step(double delay) {
		dx = ax;
		dy = ay;
		x += dx * delay;
		y += dy * delay;
		
		ax = 0;
		ay = 0;
		
	}
	
	@Override
	public void draw(Graphics g) {
		g.setColor(Color.BLACK);
		super.draw(g);
		g.setColor(Color.RED);
		for(int i = 0; i < xLoad.size(); i++) {
			double axx = xLoad.get(i);
			double ayy = yLoad.get(i);
			double accelaration = Math.sqrt(Math.pow(axx,2)+Math.pow(ayy,2));	
			drawArrowLine(g,(int)x,(int)y,(int)(x+arrowScale*axx),(int)(y+arrowScale*ayy),(int)(arrowScale*accelaration/10), (int)(arrowScale*accelaration/20));
		}
		g.setColor(Color.GREEN);
		for(int i = 0; i < xLoadGreen.size(); i++) {
			double axx = xLoadGreen.get(i);
			double ayy = yLoadGreen.get(i);
			double accelaration = Math.sqrt(Math.pow(axx,2)+Math.pow(ayy,2));	
			drawArrowLine(g,(int)x,(int)y,(int)(x+arrowScale*axx),(int)(y+arrowScale*ayy),(int)(arrowScale*accelaration/10), (int)(arrowScale*accelaration/20));
		}
		g.setColor(Color.BLUE);
		for(int i = 0; i < xLoadBlue.size(); i++) {
			double axx = xLoadBlue.get(i);
			double ayy = yLoadBlue.get(i);
			double accelaration = Math.sqrt(Math.pow(axx,2)+Math.pow(ayy,2));	
			drawArrowLine(g,(int)x,(int)y,(int)(x+arrowScale*axx),(int)(y+arrowScale*ayy),(int)(arrowScale*accelaration/10), (int)(arrowScale*accelaration/20));
		}
		g.setColor(Color.CYAN);
		for(int i = 0; i < xLoadCyan.size(); i++) {
			double axx = xLoadCyan.get(i);
			double ayy = yLoadCyan.get(i);
			double accelaration = Math.sqrt(Math.pow(axx,2)+Math.pow(ayy,2));	
			drawArrowLine(g,(int)x,(int)y,(int)(x+arrowScale*axx),(int)(y+arrowScale*ayy),(int)(arrowScale*accelaration/10), (int)(arrowScale*accelaration/20));
		}
		
		xLoad = new ArrayList<Double>();
		yLoad = new ArrayList<Double>();
		xLoadGreen = new ArrayList<Double>();
		yLoadGreen = new ArrayList<Double>();
		xLoadBlue = new ArrayList<Double>();
		yLoadBlue = new ArrayList<Double>();
		xLoadCyan = new ArrayList<Double>();
		yLoadCyan = new ArrayList<Double>();
		g.setColor(Color.BLACK);
	}
	
	/**
	 * Draw an arrow line between two points.
	 * @param g the graphics component.
	 * @param x1 x-position of first point.
	 * @param y1 y-position of first point.
	 * @param x2 x-position of second point.
	 * @param y2 y-position of second point.
	 * @param d  the width of the arrow.
	 * @param h  the height of the arrow.
	 */
	private void drawArrowLine(Graphics g, int x1, int y1, int x2, int y2, int d, int h) {
	    int dx = x2 - x1, dy = y2 - y1;
	    double D = Math.sqrt(dx*dx + dy*dy);
	    double xm = D - d, xn = xm, ym = h, yn = -h, x;
	    double sin = dy / D, cos = dx / D;
	    x = xm*cos - ym*sin + x1;
	    ym = xm*sin + ym*cos + y1;
	    xm = x;

	    x = xn*cos - yn*sin + x1;
	    yn = xn*sin + yn*cos + y1;
	    xn = x;

	    int[] xpoints = {x2, (int) xm, (int) xn};
	    int[] ypoints = {y2, (int) ym, (int) yn};

	    g.drawLine(x1, y1, x2, y2);
	    g.fillPolygon(xpoints, ypoints, 3);
	}

}

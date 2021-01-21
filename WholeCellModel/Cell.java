import java.awt.*;
import java.util.*;

public class Cell {
	public boolean hasLoad;

	{
		hasLoad = false;
	}

	// each cell is defined by an array of nodes
	private ArrayList<BlobC> blobs;
	// the center of the cell
	private double cx, cy;

	public double initialArea;

	// for the drawing of arrows to show forces
	private boolean drawElastic = false, drawDamping = false, drawOsmosis = false, drawStiff = false, drawHydrostatic = false, drawCompress = false, drawLateral = false, drawEcto = false;
	private Color dampingColor = Color.GREEN, osmosisColor = Color.blue;
	// initial length of the edges
	private ArrayList<Double> elasticLengths;

	// constructor: this function is called first when we initiate a cell object
	public Cell(ArrayList<BlobC> blobs) {
		this.blobs = blobs;
		elasticLengths = new ArrayList<Double>();
		BlobC blob1;
		BlobC blob2 = blobs.get(0);
		double x1, x2, y1, y2;
		x2 = blob2.getX();
		y2 = blob2.getY();
		for(int i = 0; i < blobs.size(); i++) {
			x1 = x2;
			y1 = y2;
			blob1 = blob2;
			blob2 = blobs.get((i+1)%blobs.size());
			x2 = blob2.getX();
			y2 = blob2.getY();
			double elasticLength = Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
			elasticLengths.add(elasticLength);
			
			initialArea = area();
			double sumX = 0, sumY = 0;
			for(BlobC blob: blobs) {
				sumX += blob.getX();
				sumY += blob.getY();
			}
			cx = sumX/blobs.size();
			cy = sumY/blobs.size();
		}
	}
	// return the ith node
	public BlobC getBlob(int i) {
		return blobs.get(i);
	}
	// find the size of the array of nodes
	public int getSize() {
		return blobs.size();
	}

	// calculate the area of the cell
	public double area() {
		double x1 = 0.0, y1 = 0.0;
		double x2 = blobs.get(0).getX();
		double y2 = blobs.get(0).getY();
		double x3 = blobs.get(1).getX();
		double y3 = blobs.get(1).getY();
		int size = blobs.size();
		double area = 0.0;
		for(int i = 0; i < size; i++) {
			x1 = x2;
			y1 = y2;
			x2 = x3;
			y2 = y3;
			x3 = blobs.get((i+2)%size).getX();
			y3 = blobs.get((i+2)%size).getY();
			area += (x2*(y1-y3));
		}
		
		return 0.5*Math.abs(area);
	}

	// apply the elastic force
	public void elastic(double elasticConstant) {
		int size = blobs.size();
		double x1 = 0.0, y1 = 0.0;
		double x2 = blobs.get(0).getX();
		double y2 = blobs.get(0).getY();
		for(int i = 0; i < size; i++) {
			x1 = x2;
			y1 = y2;
			x2 = blobs.get((i+1)%size).getX();
			y2 = blobs.get((i+1)%size).getY();
			double length = Math.sqrt(Math.pow((x1-x2), 2) + Math.pow((y1-y2),2));
			double f = elasticConstant*(length-elasticLengths.get(i));
			double fx = f*(x2-x1)/length;
			double fy = f*(y2-y1)/length;
			blobs.get(i).applyLoad(fx, fy, drawElastic);
			blobs.get((i+1)%size).applyLoad(-fx, -fy, drawElastic);
		}
//		System.out.println(elasticLengths);
	}

	// apply the volume conservation pressure
	public void osmosis(double osmosisConstant) {
		double f = osmosisConstant*(area()-initialArea);
		double[][] dirs = normal();
		for(int i = 0; i< blobs.size(); i++) {
			double[] normalDirection = dirs[i];
			blobs.get(i).applyLoad(-f*normalDirection[0], -f*normalDirection[1], drawOsmosis);
		}
	}

	// apply the yolk hydrostatic pressure
	public void hydrostatic(double hydroConstant) {
		double f = hydroConstant;
		double[][] dirs = normal();
		for(int i = 0; i < blobs.size(); i++) {
			double[] normalDirection = dirs[i];
			blobs.get(i).applyLoad(f*normalDirection[0], f*normalDirection[1], drawHydrostatic);
		}
	}

	// apply the apical constriction
	public void constrictLoad(int i, double compressConstant, double compressRatio) {
		BlobC blob1 = blobs.get((i+blobs.size())%blobs.size());
		BlobC blob2 = blobs.get((i+1+blobs.size())%blobs.size());
		double x1 = blob1.getX();
		double y1 = blob1.getY();
		double x2 = blob2.getX();
		double y2 = blob2.getY();
		double length = Math.sqrt(Math.pow((x1-x2), 2) + Math.pow((y1-y2),2));
		double f = compressConstant*(length-(compressRatio* elasticLengths.get(i)));
		double fx = f*(x2-x1)/length;
		double fy = f*(y2-y1)/length;
		blob1.applyLoad(fx, fy, drawCompress);
		blob2.applyLoad(-fx, -fy, drawCompress);
	}

	// apply the lateral constriction
	public void lateralLoad(int i, double lateralConstant, double lateralRatio) {
		BlobC blob1 = blobs.get((i+blobs.size())%blobs.size());
		BlobC blob2 = blobs.get((i+1+blobs.size())%blobs.size());
		double x1 = blob1.getX();
		double y1 = blob1.getY();
		double x2 = blob2.getX();
		double y2 = blob2.getY();
		double length = Math.sqrt(Math.pow((x1-x2), 2) + Math.pow((y1-y2),2));
		double f = lateralConstant*(length-(lateralRatio* elasticLengths.get(i)));
		double fx = f*(x2-x1)/length;
		double fy = f*(y2-y1)/length;
		blob1.applyLoad(fx, fy, drawLateral);
		blob2.applyLoad(-fx, -fy, drawLateral);
	}

	// not used
	// apply the apical-basal elongation on ectodermal cells
	public void ectoLoad(int i, double ectoConstant, double ectoRatio) {
		BlobC blob1 = blobs.get((i+blobs.size())%blobs.size());
		BlobC blob2 = blobs.get((i+1+blobs.size())%blobs.size());
		double x1 = blob1.getX();
		double y1 = blob1.getY();
		double x2 = blob2.getX();
		double y2 = blob2.getY();
		double length = Math.sqrt(Math.pow((x1-x2), 2) + Math.pow((y1-y2),2));
		double f = ectoConstant*(length-(ectoRatio* elasticLengths.get(i)));
		double fx = f*(x2-x1)/length;
		double fy = f*(y2-y1)/length;
		blob1.applyLoad(fx, fy, drawEcto);
		blob2.applyLoad(-fx, -fy, drawEcto);
	}

	// the corner stiffness load
	// to keep the corners approximately 90 degrees
	public void cornerstiff(double stiffnessConstant, int nNode11, int nNode12, int nNode2) {
		
		int size = nNode11 + nNode12 + 2*nNode2;
	//	get the three nodes on the corner
		BlobC blob1 = blobs.get(size-1);
		BlobC blob2 = blobs.get(0);
		BlobC blob3 = blobs.get(1);
		// calculate the force
		double[] force1 = stiffnessPerpendHelper(blob1.getX(),blob1.getY(),blob2.getX(),blob2.getY(),blob3.getX(),blob3.getY(),stiffnessConstant);
		blob2.applyLoad(force1[0], force1[1], drawStiff);
		
		blob1 = blobs.get(nNode11-1);
		blob2 = blobs.get(nNode11);
		blob3 = blobs.get(nNode11+1);
		double[] force2 = stiffnessPerpendHelper(blob1.getX(),blob1.getY(),blob2.getX(),blob2.getY(),blob3.getX(),blob3.getY(),stiffnessConstant);
		blob2.applyLoad(force2[0], force2[1], drawStiff);
		
		blob1 = blobs.get(nNode11+nNode2-1);
		blob2 = blobs.get(nNode11+nNode2);
		blob3 = blobs.get(nNode11+nNode2+1);
		double[] force3 = stiffnessPerpendHelper(blob1.getX(),blob1.getY(),blob2.getX(),blob2.getY(),blob3.getX(),blob3.getY(),stiffnessConstant);
		blob2.applyLoad(force3[0], force3[1], drawStiff);
		
		blob1 = blobs.get(nNode11+nNode2+nNode12-1);
		blob2 = blobs.get(nNode11+nNode2+nNode12);
		blob3 = blobs.get((nNode11+nNode2+nNode12+1)%size);
		double[] force4 = stiffnessPerpendHelper(blob1.getX(),blob1.getY(),blob2.getX(),blob2.getY(),blob3.getX(),blob3.getY(),stiffnessConstant);
		blob2.applyLoad(force4[0], force4[1], drawStiff);
	}
	
	public double[] stiffnessPerpendHelper (double x1, double y1, double x2, double y2, double x3, double y3, double k) {
		double [] res = new double[2];
		// to tell if the current angle is 0-180, or 180-360, we need to see if the three nodes are clockwise or counter-clockwise
		double clockwise = (x1*(y2-y3) + y1*(x3-x2) + x2*y3-x3*y2);
		double length = Math.sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
		double cosTheta = ((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2))/(Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)))/(Math.sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)));
		double theta;
		if(clockwise>0) {
			theta = (3*(Math.PI/2)-Math.acos(cosTheta));
			
			if((x1-x2)*(y1-y3)+(y1-y2)*(x3-x1)>0) {

				res[0] = k*theta*(y1-y3)/length;
				res[1] = k*theta*(x3-x1)/length;
			}else {

				res[0] = -k*theta*(y1-y3)/length;
				res[1] = -k*theta*(x3-x1)/length;
			}
		}
		else {

			theta = Math.acos(cosTheta)-Math.PI/2;
			if((x1-x2)*(y1-y3)+(y1-y2)*(x3-x1)<0) {
				res[0] = k*theta*(y1-y3)/length;
				res[1] = k*theta*(x3-x1)/length;
			}else {
				res[0] = -k*theta*(y1-y3)/length;
				res[1] = -k*theta*(x3-x1)/length;
			}
		}
		
		return res;
	}

	public void draw(Graphics g) {
		int size = blobs.size();
		double x1 = 0.0, y1 = 0.0;
		double x2 = blobs.get(0).getX();
		double y2 = blobs.get(0).getY();
		for(int i = 0; i < size; i++) {
			blobs.get(i).draw(g);
			x1 = x2;
			y1 = y2;
			x2 = blobs.get((i+1)%size).getX();
			y2 = blobs.get((i+1)%size).getY();
			g.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
		}
		
	}

	// to calculate the normal direction of each node of the cell
	// the normal direction points outward
	public double[][] normal() {
		double currentArea = area();
		
		int size = blobs.size();
		double[][] res = new double[size][2];
		for(int i = 0; i < size; i++) {
			double x1 = blobs.get((i-1+size)%size).getX();
			double y1 = blobs.get((i-1+size)%size).getY();
			double x2 = blobs.get((i+1)%size).getX();
			double y2 = blobs.get((i+1)%size).getY();
			double length = Math.sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
			double[] dir = {(y2-y1)/length,-(x2-x1)/length};
			Blob blob0 = blobs.get(i);
			
			
			double x3 = blob0.getX();
			double y3 = blob0.getY();
			
			blob0.setX(blob0.getX()+0.5*dir[0]);
			blob0.setY(blob0.getY()+0.5*dir[1]);
			double newArea = area();
			if(newArea>currentArea) {
				res[i][0]=dir[0];
				res[i][1]=dir[1];
			}
			else {
				res[i][0]=-dir[0];
				res[i][1]=-dir[1];
			}
			blob0.setX(blob0.getX()-0.5*dir[0]);
			blob0.setY(blob0.getY()-0.5*dir[1]);
		}
		return res;
	}

}

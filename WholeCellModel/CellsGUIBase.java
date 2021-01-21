
import java.awt.*;
import java.io.*;
import java.io.FileNotFoundException;

import javax.swing.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.lang.Math;
import java.util.*;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.time.Instant;
/**
 * Animated cells.  
 * 
 * @author Chris Bailey-Kellogg, Dartmouth CS 10, Spring 2016, based on animated agents from previous terms
 * @author Shicheng Huang
 */
public class CellsGUIBase extends DrawingGUI {
	
	protected static int width=800, height=700;	// size of the world
	protected static double r1 = 300.0, r2 = 280.0, r3 = 301.0, r4 = 340.0, scaleMultiplier = 0.95;
	
	private int nMesh = 320;
	private int nNode11 = 1, nNode12 = 1, nNode2 = 1, nNode3 = 1; // number of nodes on each edge
	
	protected static long timeCount = 0;
	protected static final int cutoffTime = 3000, stopWriting = 600000, cutoffForce = 10;
	protected int delay = 50;							// for the timer
	
	protected double elasticConstant = 0.2;
//	protected double dampingConstant = 0.1;
	protected double osmosisConstant = 0.1;
	protected double tugerOsmosisConstant = 0.001;
	protected double tugerHydroConstant = 0.05;
	protected double leonardJonesConstant = 0.001, leonardJonesR0 = 2, rMax = 50, rTruncate = 0.33;
//	protected double glueConstant = 5.0, glueR0 = 1, glueRTruncate = 0.25;
	
	protected boolean drawLeonard = false;
	protected boolean drawGlue = false;
	protected boolean drawStiff = false;
	
	private int compressNum = 36;   // the number of constricting cells
	double compressConstant[] = {0.75,1.96}, compressRatio[] = {0.01,0.01}, compressConstPara;
	double[] compressDistribution;
//	double compressConstant[] = {0.5,0.0}, compressRatio[] = {0.01,0.01};
	
	private int ectoBegin = 36, ectoEnd = 160;	// the number of lengthening ectoderm cells
	double  ectoConstant[] = {0.2,0.2,0.2,0.2}, ectoRatio[] = {1.5,1.5,1.5,1.5};
	
	double lateralConstant[] = {0.2,0.2}, lateralRatio[] = {0.01,0.01}, lateralPara[];
	
	double stiffnessConstant = 0.8;
	
	protected ArrayList<Cell> intraCells;				// list of boundary cells
	protected ArrayList<Cell> cells;					// list of all the cells to handle
	protected ArrayList<BlobC> blobCells;			// list of all the blobs in epithelial cells
	protected Cell tugerCell;							// the inner tuger cell
	
	BufferedWriter writer;
	
	public CellsGUIBase() {
		super("Animated Cells", width, height);
		blobCells = new ArrayList<BlobC>();
		cells = new ArrayList<Cell>();
		intraCells = new ArrayList<Cell>();
		
		compressDistribution = new double[compressNum*nNode11];
		lateralPara = new double[2];
		Scanner scan;
		
		String s = "cell/distribution-"+(nNode11*compressNum*2)+".txt";
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
		String ts  = dateFormat.format(timestamp);
		String s2 = "Gaussian-"+compressNum+"-stiff-"+stiffnessConstant+"-lateral-"+lateralConstant[0]+"-"+lateralRatio[0]+"-hydro-"+tugerHydroConstant+"-osmosis-"+osmosisConstant+"-cutoffTime-"+cutoffTime+"-cutoffForce-"+cutoffForce+"-"+ts+".txt";
//		String s2 = "Gaussian-"+compressNum+"-stiff-"+stiffnessConstant+"-lateral-"+lateralConstant[0]+'-'+lateralRatio[0]+"-hydro-"+tugerHydroConstant+"-osmosis-"+osmosisConstant+"-cutoffTime-"+cutoffTime+"-cutoffForce-"+cutoffForce+".txt";
		try {
			File file = new File(s);
			scan = new Scanner(file);
			writer = new BufferedWriter(new FileWriter(s2,false));
			for(int i = 0; i < compressNum*nNode11; i++) {
				compressDistribution[i] = scan.nextDouble();
				System.out.println(compressDistribution[i]);
			}
			scan.close();
			String s3 = "time        node location";
			writer.write(s3);
			writer.newLine();
		}
		catch(IOException e1){
			e1.printStackTrace();
		}
		
		try {
			meshing();
//			scale(scaleMultiplier);
			
		}
		catch (Exception e){
			System.out.println("error in meshing!");
		}
		
		
		// Timer drives the animation.
		startTimer();
	}
//	
	
	public void meshing() {
		List<List<BlobC>> lateral = new ArrayList<>();
		List<List<BlobC>> intraLateral = new ArrayList<>();
		List<BlobC> basel = new ArrayList<>();
		List<BlobC> apical = new ArrayList<>();
		List<BlobC> intraApical = new ArrayList<>();
		List<BlobC> intraBasel = new ArrayList<>();
		double s1, c1;
		
		for(int i = 0; i < nMesh; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh));
			List<BlobC> temp = new ArrayList<>();
			for(int j = 1; j < nNode2; j++) {
				double r = r1 + (r2-r1)/nNode2*j;
				BlobC bTemp = new BlobC(r*s1+width/2,r*c1+height/2);
				temp.add(bTemp);
			}
			lateral.add(temp);
			blobCells.addAll(temp);
			temp = new ArrayList<>();
			for(int j = 1; j < nNode3; j++) {
				double r = r4 + (r3-r4)/nNode3*j;
				BlobC bTemp = new BlobC(r*s1+width/2,r*c1+height/2);
				temp.add(bTemp);
			}
			intraLateral.add(temp);
		}
		for(int i = 0; i < nMesh*nNode11; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode11));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode11));
			apical.add(new BlobC(r1*s1+width/2,r1*c1+height/2));
		}
		blobCells.addAll(apical);
		for(int i = 0; i < nMesh*nNode12; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode12));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode12));
			basel.add(new BlobC(r2*s1+width/2,r2*c1+height/2));
		}
		blobCells.addAll(basel);
		for(int i = 0; i < nMesh*nNode12; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode12));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode12));
			intraApical.add(new BlobC(r4*s1+width/2,r4*c1+height/2));
		}
		for(int i = 0; i < nMesh*nNode12; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode12));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode12));
			intraBasel.add(new BlobC(r3*s1+width/2,r3*c1+height/2));
		}
		
		for(int i = 0; i < nMesh; i++) {
			ArrayList<BlobC> blobs = new ArrayList<BlobC>();
			for(int j = 0; j < nNode11; j++) {
				BlobC blob = apical.get(i*nNode11+j);
				blobs.add(blob);
			}
			
			blobs.add(apical.get(((i+1)*nNode11)%(nNode11*nMesh)));
			
			for(int j = 0; j < nNode2-1; j++) {
				BlobC blob = lateral.get((i+1)%nMesh).get(j);
				blobs.add(blob);
			}
			
			blobs.add(basel.get(((i+1)*nNode12)%(nNode12*nMesh)));
			
			for(int j = 0; j < nNode12; j++) {
				BlobC blob = basel.get((i+1)*nNode12-j-1);
				blobs.add(blob);
			}
			
			for(int j = 0; j < nNode2-1; j++) {
				BlobC blob = lateral.get(i).get(nNode2-j-2);
				blobs.add(blob);
			}
			
			Cell cell = new Cell(blobs);
			cells.add(cell);
		}
		
		for(int i = 0; i < nMesh; i++) {
			ArrayList<BlobC> blobs = new ArrayList<>();
			blobs.addAll(intraApical.subList(i*nNode12,(i+1)*nNode12));
			blobs.add(intraApical.get(((i+1)*nNode12)%(nNode12*nMesh)));
			blobs.addAll(intraLateral.get((i+1)%nMesh));
			
			List<BlobC> tempBasel = intraBasel.subList(i*nNode12,(i+1)*nNode12);
			
			blobs.add(intraBasel.get(((i+1)*nNode12)%(nNode12*nMesh)));
			Collections.reverse(tempBasel);
			blobs.addAll(tempBasel);
			Collections.reverse(tempBasel);
			List<BlobC> tempLateral = intraLateral.get(i);
			Collections.reverse(tempLateral);
			blobs.addAll(tempLateral);
			Collections.reverse(tempLateral);
			Cell cell = new Cell(blobs);
			intraCells.add(cell);
		}
//		System.out.println("start tuger!");
		tugerCell = new Cell((ArrayList<BlobC>)basel);
//		tugerCell.initialArea = 1.05*tugerCell.initialArea;
				
	}
	
//	public void scale(double multiplier) {
//		for(Cell cell: cells) {
//			cell.scale(multiplier);
//		}
////		tugerCell.scale(multiplier);
//	}
	

	@Override
	public void handleKeyPress(char k) {
		System.out.println("Handling key '"+k+"'");
		if (k == 'f') { // faster
			if (delay>1) delay /= 2;
			setTimerDelay(delay);
			System.out.println("delay:"+delay);
		}
		else if (k == 's') { // slower
			delay *= 2;
			setTimerDelay(delay);
			System.out.println("delay:"+delay);
		}
		else { // blob type
			// handle measurement		
		}
	}
	
	/**
	 * DrawingGUI method, here just drawing all the cells
	 */
	@Override
	public void draw(Graphics g) {
		// Ask all the cells to draw themselves.
		for (Cell cell : cells) {
			cell.draw(g);
		}		
		for (Cell cell: intraCells) {
			cell.draw(g);
		}
		tugerCell.draw(g);
	}
	
	public void timeDependent() {
		timeCount++;
		if(timeCount <= cutoffTime) {
			double frac = (cutoffTime-timeCount)*1.0/cutoffTime;
			double factor1 = cutoffForce - (cutoffForce-1)*frac*frac;
			double factor21 = lateralRatio[0] + (1.0-lateralRatio[0])*frac;
			double factor22 = lateralRatio[1] + (1.0-lateralRatio[1])*frac;
			System.out.println("Apical constrction-"+factor1);
			compressConstPara = compressConstant[0]*factor1;
			
			lateralPara[0] = factor21;
			lateralPara[1] = factor22;
			
//			compressConstant[0] = 1.001 * compressConstant[0];
//			compressConstant[0] = 10;
//			System.out.println("Apical constrction-"+compressConstant[0]);
		}	
	}
	
	public void leonardJones() {
		PointQuadtree<BlobC> blobTree = new PointQuadtree<BlobC>(blobCells.get(0),0,0,width,height);
		for(int i = 1; i < blobCells.size(); i++) {
			blobTree.insert(blobCells.get(i));
		}
		
		for(int j = 0; j < intraCells.size(); j++) {
			
			Cell intraCell = intraCells.get(j);
			BlobC blob1, blob2;
			blob2 = intraCell.getBlob(0);
			for(int i = 0; i < intraCell.getSize(); i++) {
				blob1 = blob2;
				blob2 = intraCell.getBlob((i+1)%intraCell.getSize());
				List<BlobC> blobTrack = new ArrayList<>();
				blobTrack = blobTree.findInCircle((blob1.getX()+blob2.getX())/2, (blob1.getY()+blob2.getY())/2, rMax);
				for(BlobC blob : blobTrack) {
					if(blob == blob1 || blob == blob2) continue;
					double f[] = leonardJonesHelper(blob1,blob2,blob,leonardJonesConstant,leonardJonesR0);
					blob.applyLoad(f[0], f[1], drawLeonard);
				}
			}
		}
		
		for(int j = 0; j < cells.size(); j++) {
			Cell cell = cells.get(j);
			BlobC blob1, blob2;
			blob2 = cell.getBlob(0);
			for(int i = 0; i < nNode11 + nNode2 + nNode12; i++) {
				blob1 = blob2;
				blob2 = cell.getBlob(i+1);
				List<BlobC> blobTrack = new ArrayList<>();
				blobTrack = blobTree.findInCircle((blob1.getX()+blob2.getX())/2, (blob1.getY()+blob2.getY())/2, rMax);
				for(BlobC blob : blobTrack) {
					if(blob == blob1 || blob == blob2) continue;
					double f[] = leonardJonesHelper(blob1,blob2,blob,leonardJonesConstant,leonardJonesR0);
					blob.applyLoad(f[0], f[1], drawLeonard);
				}
			}
		}
	}
	

	
	public double[] leonardJonesHelper(Blob blob1, Blob blob2, Blob blob3, double kij, double r0) {
		double x1 = blob1.getX(), y1 = blob1.getY();
		double x2 = blob2.getX(), y2 = blob2.getY();
		double x3 = blob3.getX(), y3 = blob3.getY();
		double rij = 0.0;
		double length = Math.sqrt(Math.pow(x2-x1,2)+Math.pow(y2-y1,2));
		double proj1 = (x3-x1)*(x2-x1)/length+(y3-y1)*(y2-y1)/length;
		double proj2 = (x3-x2)*(x2-x1)/length+(y3-y2)*(y2-y1)/length;
		double f = 0.0;
		double[] direction = {0.0, 0.0};
		if(Math.abs(proj1) > Math.abs(proj2)) {
			if (proj2 > 0) {
				rij = Math.sqrt((Math.pow(x3-x2, 2)+Math.pow(y3-y2, 2)));
				direction[0] = (x3-x2)/rij; direction[1] = (y3-y2)/rij;
			}else {
				rij = Math.sqrt((Math.pow(x3-x1, 2)+Math.pow(y3-y1, 2)-Math.pow(proj1, 2)));
				if((x3-x1)*(y1-y2)+(y3-y1)*(x2-x1)>0){
					direction[0]=(y1-y2)/length; direction[1] = (x2-x1)/length;
				}else {
					direction[0]=(y2-y1)/length; direction[1] = (x1-x2)/length;
				}
			}
		}else {
			if(proj1 > 0) {
				rij = Math.sqrt((Math.pow(x3-x1, 2)+Math.pow(y3-y1, 2)-Math.pow(proj1, 2)));
				if((x3-x1)*(y1-y2)+(y3-y1)*(x2-x1)>0){
					direction[0]=(y1-y2)/length; direction[1] = (x2-x1)/length;
				}else {
					direction[0]=(y2-y1)/length; direction[1] = (x1-x2)/length;
				}
			}else {
				rij = Math.sqrt((Math.pow(x3-x1, 2)+Math.pow(y3-y1, 2)));
				direction[0] = (x3-x1)/rij; direction[1] = (y3-y1)/rij;
			}
		}
		if(rij < rTruncate) rij = rTruncate;
//		f = kij/Math.pow(rij, 3)*(8*Math.pow(r0/rij, 6)-10*Math.pow(r0/rij, 3));
//		System.out.println("rij "+rij);
		f = kij/Math.pow(rij, 3)*(6*Math.pow(r0/rij, 4)-8*Math.pow(r0/rij, 2));
//		System.out.println(f);
//		System.out.println("leonard "+f);
		double ff[] = {f*direction[0], f*direction[1]};
		return ff;
	}
	
//	public void cornerStiff() {
//		Cell cell2 = cells.get(0);
//		for(int i = 0; i <compressNum+1; i++) {
//			Cell cell1 = cell2;
//			cell2 = cells.get(i+1);
//			BlobCell blob11 = cell1.getBlob(nNode11+nNode2+1);
//			BlobCell blob12 = cell1.getBlob(nNode11+nNode2);
//			BlobCell blob21 = cell2.getBlob(nNode11+nNode12+nNode2);
//			BlobCell blob22 = cell2.getBlob(nNode11+nNode12+nNode2-1);
//			double[] force = stiffnessHelper(blob11.getX(),blob11.getY(),blob12.getX(),blob12.getY(),blob22.getX(),blob22.getY(),stiffnessConstant);
//			blob12.applyLoad(force[0], force[1], drawStiff);
//			force = stiffnessHelper(blob11.getX(),blob11.getY(),blob21.getX(),blob21.getY(),blob22.getX(),blob22.getY(),stiffnessConstant);
//			blob21.applyLoad(force[0], force[1], drawStiff);
//			
//		}
//		cell2 = cells.get(cells.size()-compressNum-1-1);
//		for(int i = cells.size()-compressNum-1-1; i <cells.size(); i++) {
//			Cell cell1 = cell2;
//			cell2 = cells.get((i+1)%cells.size());
//			BlobCell blob11 = cell1.getBlob(nNode11+nNode2+1);
//			BlobCell blob12 = cell1.getBlob(nNode11+nNode2);
//			BlobCell blob21 = cell2.getBlob(nNode11+nNode12+nNode2);
//			BlobCell blob22 = cell2.getBlob(nNode11+nNode12+nNode2-1);
//			double[] force = stiffnessHelper(blob11.getX(),blob11.getY(),blob12.getX(),blob12.getY(),blob22.getX(),blob22.getY(),stiffnessConstant);
//			blob12.applyLoad(force[0], force[1], drawStiff);
//			force = stiffnessHelper(blob11.getX(),blob11.getY(),blob21.getX(),blob21.getY(),blob22.getX(),blob22.getY(),stiffnessConstant);
//			blob21.applyLoad(force[0], force[1], drawStiff);
//			
//		}
//	}
	
	public double[] stiffnessHelper(double x1, double y1, double x2, double y2, double x3, double y3, double k) {
		double length = Math.abs(((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1))/(Math.sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1))));
		double normal = Math.sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1));
		if((x2-x1)*(y1-y3)+(y2-y1)*(x3-x1)>0) {
			double[] res = {-k*length*(y1-y3)/normal,-k*length*(x3-x1)/normal};
			return res;
		}
		double[] res = {k*length*(y1-y3)/normal,k*length*(x3-x1)/normal};
		return res;
	}
	
	public void ectoLoad() {
		double i1 = ectoRatio[0];
		double i2 = ectoRatio[1];
		double i3 = ectoRatio[2];
		double i4 = ectoRatio[3];
		double j1 = ectoConstant[0];
		double j2 = ectoConstant[1];
		double j3 = ectoConstant[2];
		double j4 = ectoConstant[3];
		System.out.println("ecto-"+j1+"-"+i1);
		
		int ectoNum = ectoEnd-ectoBegin;
		double func1[] = new double[ectoNum*nNode11];
		double func2[] = new double[ectoNum*nNode11];
		for(int i = 0; i < ectoNum*nNode11; i++) {
			func1[i] = (i2-i1)*(1.0/ectoNum/nNode11*i)+i1;
		}
		for(int j = 0; j < ectoNum*nNode11; j++) {
			func2[j] = (j2-j1)*(1.0/ectoNum/nNode11*j)+j1;
		}
		int count = 0;
		for(int i = ectoBegin; i < ectoEnd; i++) {
			Cell cell = cells.get(i);
			for(int j = 0; j < nNode11; j++) {
				cell.ectoLoad(j,func2[count],func1[count]);
				count++;
			}
		}
		count = 0;
		for(int i = nMesh - ectoBegin-1; i >= nMesh - ectoEnd; i--) {
			Cell cell = cells.get(i);
			for(int j = nNode11-1; j >=0; j--) {
				cell.ectoLoad(j,func2[count],func1[count]);
				count++;
			}
		}
		
		func1 = new double[ectoNum*nNode12];
		func2 = new double[ectoNum*nNode12];
		for(int i = 0; i < ectoNum*nNode12; i++) {
			func1[i] = (i4-i3)*(1.0/ectoNum/nNode12*i)+i3;
		}
		for(int j = 0; j < ectoNum*nNode11; j++) {
			func2[j] = (j4-j3)*(1.0/ectoNum/nNode12*j)+j3;
		}
		count = 0;
		for(int i = ectoBegin; i < ectoEnd; i++) {
			Cell cell = cells.get(i);
			for(int j = nNode11+nNode2+nNode12-1; j >= nNode11+nNode2; j--) {
				cell.ectoLoad(j,func2[count],func1[count]);
				count++;
			}
		}
		count = 0;
		for(int i = nMesh - ectoBegin-1; i >= nMesh - ectoEnd; i--) {
			Cell cell = cells.get(i);
			for(int j = nNode11+nNode2; j < nNode11 + nNode2+ nNode12; j++) {
				cell.ectoLoad(j,func2[count],func1[count]);
				count++;
			}
		}
	}
	
	
	public void lateralLoad() {
//		double i1 = lateralPara[0];
//		double i2 = lateralPara[1];
		double i1 = lateralRatio[0];
		double i2 = lateralRatio[1];
		double j1 = lateralConstant[0];
		double j2 = lateralConstant[1];
		System.out.println("lateral-"+j1+"-"+i1);
		double func1[] = new double[nNode2];
		double func2[] = new double[nNode2];
		for(int i = 0; i < nNode2; i++) {
			func1[i] = (i2-i1)*(1.0/nNode2*i)+i1;
		}
		for(int j = 0; j < nNode2; j++) {
			func2[j] = (j2-j1)*(1.0/nNode2*j)+j1;
		}
		for(int i = 0; i <nMesh; i++) {
			Cell cell = cells.get(i);
			for(int j = nNode11; j < nNode11+nNode2; j++) {
				cell.lateralLoad(j, func2[j-nNode11], func1[j-nNode11]);
			}
//			for(int j = nNode11+nNode12+nNode2; j < nNode11+2*nNode2+nNode12; j++) {
//				cell.compressiveLoad(j,func2[j-nNode11-nNode12-nNode2], func1[j-nNode11-nNode12-nNode2]);
//			}
		}
	}
	
	public void constrictLoad() {
		double i1 = compressRatio[0];
		double i2 = compressRatio[1];
		double j1 = compressConstant[0];
		double j2 = compressConstant[1];
		
		j1 = compressConstPara;
		
		double[] func1 = new double[compressNum*nNode11];
		double[] func2 = new double[compressNum*nNode11];
		for(int i = 0; i < compressNum*nNode11; i++) {
			func1[i] = (i2-i1)*(1.0/compressNum/nNode11*i)+i1;
		}
		double cutoff = j1* Math.exp(-0.5*j2*j2*(compressNum*nNode11-1)/nNode11/compressNum);
		for(int j = 0; j < compressNum*nNode11; j++) {
//			func2[j] = j1*compressDistribution[j];
			double loc = j2/compressNum/nNode11*j;
			
			func2[j] = j1* Math.exp(-0.5*loc*loc)-cutoff;
//			func2[j] = 0.75*func2[j];
		}
		

//		func2 is compress constant
//		func1 is compress ratio
		int count = 0;
		for(int i = 0; i<compressNum; i++) {
			Cell cell = cells.get(i);
			for(int j = 0; j<nNode11; j++) {
				cell.constrictLoad(j, func2[count], func1[count]);
//				cell.compressiveLoad(j, compressConstant, compressRatio[i]);
				count++;
			}
		}
		count = 0;
		for(int i = cells.size()-1; i>=cells.size()-compressNum; i--) {
			Cell cell = cells.get(i);
			for(int j = nNode11-1; j>=0; j--) {
				cell.constrictLoad(j, func2[count], func1[count]);
//				cell.compressiveLoad(j, compressConstant, compressRatio[cells.size()-1-i]);
				count++;
			}
		}
		
	}
	/**
	 * DrawingGUI method, here having all the cells take a step
	 */
	
	public void writeFile() {
		if(timeCount % 100 == 0) {
			int size = nNode11 + 2* nNode2 +nNode12;
			StringBuilder strb = new StringBuilder();
			strb.append(timeCount);
			strb.append("      ");
			for(Cell cell : cells) {
				for(int i = 0; i <= nNode2; i++) {
					double x1 = cell.getBlob((size-i)%size).getX();
					double y1 = cell.getBlob((size-i)%size).getY();
					if(Double.isNaN(x1)) {
						try {
							writer.close();
						}catch(IOException e1) {
							e1.printStackTrace();
						}
					}
					strb.append(x1);
					strb.append(' ');
					strb.append(y1);
					strb.append(' ');
					
				}
				
			}
			String s = strb.toString();
			try {
				writer.write(s);
				writer.newLine();
			}catch(IOException e1){
				e1.printStackTrace();
			}
			
			
		}
		if(timeCount==stopWriting) {
			try {
				writer.close();
			}catch(IOException e1) {
				e1.printStackTrace();
			}
			
		}
	}
	
	@Override
	public void handleTimer() {
		timeDependent();
		// Ask all the cells to move themselves.
		for (Cell cell : cells) {
			cell.elastic(elasticConstant);		// cell elastic force (1)
//			cell.damping(dampingConstant);		// cell damping force (2)
			cell.osmosis(osmosisConstant);		// volume conservation (3)
//			cell.stiffness(stiffnessConstant,nNode1,nNode2);
			cell.cornerstiff(stiffnessConstant, nNode11, nNode12, nNode2);
			
		}
		leonardJones();	 //Leonard-Jones force (4)
		constrictLoad(); //apical construction
		lateralLoad();
//		ectoLoad();
		
//		cornerStiff();
//		tugerCell.osmosis(tugerOsmosisConstant);
		tugerCell.hydrostatic(tugerHydroConstant);
		for (BlobC blob : blobCells) {
			blob.step(delay/10000.0);
		}
		writeFile();
		// Now update the GUI.
		repaint();
	}

	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				new CellsGUIBase();
			}
		});
	}
}

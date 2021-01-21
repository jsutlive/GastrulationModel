
import begginersbook.com.WriteFileAtTimePoint;

import java.awt.*;
import java.io.*;
import java.io.FileWriter;
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
public class CellsGUI extends DrawingGUI {
	
	// width and height of the display screen
	protected static int width=800, height=700;	
	
	// r1 is the outer radius of the epithelial cells
	// r2 is the inner radius of the epithelial cells
	// r3 is the inner radius of the wall cells
	// r4 is the outer radius of the wall cells
	protected static double r1 = 300.0, r2 = 200.0, r3 = 301.0, r4 = 340.0;
	
	// nMesh is number of Mesh, which is the number of total cells
	private int nMesh = 160;
	
	// nNode11 is the number of segments on the apical edge of epithelial cells
	// nNode12 is the number of segments on the basal edge of epithelial cells and wall cells
	// nNode2 is the number of segments on the lateral edge of epithelial cells
	// nNode3 is the number of segments on the lateral edge of wall cells
	private int nNode11 = 1, nNode12 = 1, nNode2 = 4, nNode3 = 1; // number of segments on each edge

	// each time the timer ticks, timeCount increase by 1. We use timeCount to apply time-dependent forces
	protected static long timeCount = 0;
	// apicalTime is the time step for apical constriction to increase from zero to apicalForce
	// the increase of apical constriction can have different temporal profiles.
	// stopWriting is the time that we stop writing to the output txt file in case the file is too large
	// lateralTime is the time for the lateral constriction resting length to change from 1 to 0.
	protected static final int apicalTime = 3000, lateralTime = 6000, stopWriting = 600000, apicalForce = 12;

	// the time interval that the timer ticks each time, the unit is in milliseconds
	protected int delay = 50;							// for the timer

	// for cell elastic membrane forces
	protected double elasticConstant = 0.5;
	// for cell volume conservation forces
	// raising this increases odds the cells will overlap, lowering it reduces the speed in which the cells can move.
	protected double osmosisConstant = 0.015;
	// for the yolk volume conservation forces
	protected double yolkOsmosisConstant = 0.2;
	// for the yolk hydrostatic pressure
	protected double yolkHydroConstant = 0.3;
	// for the LJ forces
	// rMax is the maximum distance that we consider LJ forces between an edge and a node
	// rTruncate is the minimum distance for LJ forces. If the distance between an edge and a node is smaller than rTruncate, we use rTruncate as the distance
	protected double leonardJonesConstant = 0.001, leonardJonesR0 = 2, rMax = 50, rTruncate = 0.33;
	// for glue forces, R0 is the resting length
	// glueDist is the maximum distance that we consider glue forces between two nodes.
	//If the distance between an edge and a node is smaller than glueTruncate, we use glueTruncate as the distance
	protected double glueConstant = 3, glueR0 = 1, glueRTruncate = 0.0, glueDist = 2;

	// this is for drawing arrows to show forces
	protected boolean drawLeonard = false;
	protected boolean drawGlue = false;
	protected boolean drawStiff = false;
	protected boolean drawApicalGlue = true;

	// the number of constricting cells on each half of embryo
	private int constrictNum = 16;
	// we can use different spatial distributions. For Gaussian distribution,
	double constrictConstant[] = {0.75,1.96}, constrictRatio[] = {0.01,0.01}, constrictConstPara;
	double[] constrictDistribution;
//	double compressConstant[] = {0.5,0.0}, compressRatio[] = {0.01,0.01};
	
	private int ectoBegin = 16, ectoEnd = nMesh/2-ectoBegin;	// the number of lengthening ectoderm cells
	//ectoConstant [0,1] apical cells, [2,3] basal cells (spring constant for these values)
	double  ectoConstant[] = {0.1,0.1,0.2,0.2}, ectoRatio[] = {1,1,1,1};// originally set to 1.5
	
	double lateralConstant[] = {0.4,0.4}, lateralRatio[] = {0.01,0.01}, lateralPara[];
//	int lateralBegin = compressNum, lateralEnd = nMesh-compressNum;
	int lateralBegin = 0, lateralEnd = nMesh;

	// for the corner stiffness force
	double stiffnessConstant = 1.8;//.8 default

	// array of wall cells
	protected ArrayList<Cell> wallCells;
	// array of all the epithelial cells
	protected ArrayList<Cell> cells;
	// array of all the blobs in epithelial cells
	protected ArrayList<BlobC> blobCs;
	// the inner yolk cell
	protected Cell yolkCell;

	// tree structure to detect the blobs that are close to an edge
	protected PointQuadtree<BlobC> blobTree;
	protected PointQuadtree<BlobC> apicalTree;
	protected ArrayList<BlobC> apical;
	
	
	BufferedWriter writer;

	// constructor: when the main function starts running, this function runs first.
	public CellsGUI() {
		// call the constructor of the inherited class, in this case, DrawingGUI
		super("Animated Cells", width, height);
		// if we want to initialize variables, we need to put a key word "new" in front to tell the computer to allocate memory
		blobCs = new ArrayList<BlobC>();
		cells = new ArrayList<Cell>();
		wallCells = new ArrayList<Cell>();

		constrictDistribution = new double[constrictNum*nNode11];
		lateralPara = new double[2];
		Scanner scan;

		// read from the experimental spatial distribution of apical constriction
		String s = "distribution-"+(nNode11*constrictNum*2)+".txt";
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
		String ts  = dateFormat.format(timestamp);
		// set the name of the output file
		String s2 = "Gaussian-"+constrictNum+"-stiff-"+stiffnessConstant+"-partialLateral-"+lateralConstant[0]+"-"+lateralRatio[0]+"-hydro-"+yolkHydroConstant+"-yolkOsmosis"+yolkOsmosisConstant+"-elastic-"+elasticConstant+"-osmosis-"+osmosisConstant+"-LJ-"+leonardJonesConstant+"-"+leonardJonesR0+"-"+rMax+"-"+rTruncate+"-apicalTime-"+apicalTime+"-lateralTime-"+lateralTime+"-totalForce-"+apicalForce+"-"+ts+".txt";
//		String s2 = "Gaussian-"+compressNum+"-stiff-"+stiffnessConstant+"-lateral-"+lateralConstant[0]+'-'+lateralRatio[0]+"-hydro-"+tugerHydroConstant+"-osmosis-"+osmosisConstant+"-cutoffTime-"+cutoffTime+"-cutoffForce-"+cutoffForce+".txt";
		try {
			File file = new File(s);
			// read from spatial distribution
			scan = new Scanner(file);
			// write to output
			writer = new BufferedWriter(new FileWriter(s2,false));
			for(int i = 0; i < constrictNum*nNode11; i++) {
				constrictDistribution[i] = scan.nextDouble();
				System.out.println(constrictDistribution[i]);
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
			// generate mesh
			meshing();
//			scale(scaleMultiplier);
			
		}
		catch (Exception e){
			System.out.println("error in meshing!");
		}
		
		
		// Timer drives the animation.
		startTimer();
	}

	// generate mesh
	public void meshing() {
		// 2D array of lateral nodes
		// each subarray represents a lateral membrane
		List<List<BlobC>> lateral = new ArrayList<>();
		List<List<BlobC>> wallLateral = new ArrayList<>();
		// array of basal nodes
		List<BlobC> basel = new ArrayList<>();
		// array of apical nodes
		apical = new ArrayList<>();

		List<BlobC> wallApical = new ArrayList<>();
		List<BlobC> wallBasel = new ArrayList<>();
		double s1, c1;
		
		
		// generate mesh for lateral membrane
		for(int i = 0; i < nMesh; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh));
			List<BlobC> temp = new ArrayList<>();
			for(int j = 1; j < nNode2; j++) {
				// temp is the subarray of each lateral membrane
				double r = r1 + (r2-r1)/nNode2*j;
				BlobC bTemp = new BlobC(r*s1+width/2,r*c1+height/2);
				temp.add(bTemp);
			}
			lateral.add(temp);

			blobCs.addAll(temp);
			temp = new ArrayList<BlobC>();
			for(int j = 1; j < nNode3; j++) {
				double r = r4 + (r3-r4)/nNode3*j;
				BlobC bTemp = new BlobC(r*s1+width/2,r*c1+height/2);
				temp.add(bTemp);
			}
			wallLateral.add(temp);
		}

		// generate mesh for apical membrane
		for(int i = 0; i < nMesh*nNode11; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode11));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode11));
			apical.add(new BlobC(r1*s1+width/2,r1*c1+height/2));
		}
		blobCs.addAll(apical);
		
		// generate mesh for basal membrane
		for(int i = 0; i < nMesh*nNode12; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode12));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode12));
			basel.add(new BlobC(r2*s1+width/2,r2*c1+height/2));
		}
		blobCs.addAll(basel);

		// for wall cells, similar
		for(int i = 0; i < nMesh*nNode12; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode12));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode12));
			wallApical.add(new BlobC(r4*s1+width/2,r4*c1+height/2));
		}
		for(int i = 0; i < nMesh*nNode12; i++) {
			s1 = Math.sin((double)(i*2.0*Math.PI/nMesh/nNode12));
			c1 = Math.cos((double)(i*2.0*Math.PI/nMesh/nNode12));
			wallBasel.add(new BlobC(r3*s1+width/2,r3*c1+height/2));
		}

		// assemble the apical, lateral, basal membranes to form a cell
		for(int i = 0; i < nMesh; i++) {
			// initiate an array
			ArrayList<BlobC> blobs = new ArrayList<>();
			// put apical nodes in the array
			for(int j = 0; j < nNode11; j++) {
				BlobC blob = apical.get(i*nNode11+j);
				blobs.add(blob);
			}
			
			blobs.add(apical.get(((i+1)*nNode11)%(nNode11*nMesh)));
			// put lateral nodes in the array
			for(int j = 0; j < nNode2-1; j++) {
				BlobC blob = lateral.get((i+1)%nMesh).get(j);
				blobs.add(blob);
			}
			// put basal nodes in the array
			blobs.add(basel.get(((i+1)*nNode12)%(nNode12*nMesh)));

			for(int j = 0; j < nNode12; j++) {
				BlobC blob = basel.get((i+1)*nNode12-j-1);
				blobs.add(blob);
			}
			// put lateral nodes in the array
			for(int j = 0; j < nNode2-1; j++) {
				BlobC blob = lateral.get(i).get(nNode2-j-2);
				blobs.add(blob);
			}
			// initiate cell according to the array of blobs
			Cell cell = new Cell(blobs);
			// put the new cell in the cell array
			cells.add(cell);
		}

		// assemble wall apical, basal, lateral membrane to form a wall cell
		for(int i = 0; i < nMesh; i++) {
			ArrayList<BlobC> blobs = new ArrayList<>();
			blobs.addAll(wallApical.subList(i*nNode12,(i+1)*nNode12));
			blobs.add(wallApical.get(((i+1)*nNode12)%(nNode12*nMesh)));
			blobs.addAll(wallLateral.get((i+1)%nMesh));
			
			List<BlobC> tempBasel = wallBasel.subList(i*nNode12,(i+1)*nNode12);
			
			blobs.add(wallBasel.get(((i+1)*nNode12)%(nNode12*nMesh)));
			Collections.reverse(tempBasel);
			blobs.addAll(tempBasel);
			Collections.reverse(tempBasel);
			List<BlobC> tempLateral = wallLateral.get(i);
			Collections.reverse(tempLateral);
			blobs.addAll(tempLateral);
			Collections.reverse(tempLateral);
			Cell cell = new Cell(blobs);
			wallCells.add(cell);
		}

		yolkCell = new Cell((ArrayList<BlobC>)basel);

		
				
	}
	

	
	// not used
	// we can change the time interval of the timer by press f(faster) or s(slower) button
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
		for (Cell cell: wallCells) {
			cell.draw(g);
		}
	}

	// this function handles the temporal profile of loading
	public void timeDependent() {
		timeCount++;
		if(timeCount <= apicalTime) {
//			double frac = (cutoffTime-timeCount)*1.0/cutoffTime;
//			double factor1 = cutoffForce - (cutoffForce-1)*frac*frac;
//			double factor21 = lateralRatio[0] + (1.0-lateralRatio[0])*frac*frac;
//			double factor22 = lateralRatio[1] + (1.0-lateralRatio[1])*frac*frac;
//			System.out.println("Apical constrction-"+factor1);
//			compressConstPara = compressConstant[0]*factor1;
//			
//			lateralPara[0] = factor21;
//			lateralPara[1] = factor22;


			// calculate the apical constriction temporal profile
			double frac = timeCount*1.0/apicalTime;
			double factor1 = 1.0 + (apicalForce-1.0)*frac*frac*frac*(10.0-15.0*frac+6.0*frac*frac);
			constrictConstPara = constrictConstant[0]*factor1;
			
//			compressConstant[0] = 1.001 * compressConstant[0];
//			compressConstant[0] = 10;
//			System.out.println("Apical constrction-"+compressConstant[0]);
		}	

		// calculates the lateral constriction temporal profile
		if(timeCount <= lateralTime) {
			double frac = timeCount*1.0/lateralTime;
			lateralPara[0] = 1.0 + (lateralRatio[0]-1.0)*frac*frac*frac*(10.0-15.0*frac+6.0*frac*frac);
			lateralPara[1] = 1.0 + (lateralRatio[1]-1.0)*frac*frac*frac*(10.0-15.0*frac+6.0*frac*frac);
		}
	}

	// not used
	// to glue apical membranes if they are too close
	public void apicalGlue() {
		apicalTree = new PointQuadtree<BlobC>(apical.get(0),0,0,width,height);
		for(int i = 1; i < apical.size(); i++) {
			apicalTree.insert(apical.get(i));
		}
		
		for(BlobC blobC: apical) {
			List<BlobC> blobTrack = apicalTree.findInCircle(blobC.getX(), blobC.getY(), glueDist);
			for(BlobC blobC2: blobTrack) {
				if(blobC2 != blobC) {
					double[] f = apicalGlueHelper(blobC, blobC2);
					blobC.applyLoad(f[0], f[1], drawApicalGlue);
				}
			}
		}
	}
	// calculate the glue force given the
	public double[] apicalGlueHelper(BlobC blobCell1, BlobC blobCell2) {
		double x1 = blobCell1.getX();
		double y1 = blobCell1.getY();
		double x2 = blobCell2.getX();
		double y2 = blobCell2.getY();
		double length = Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
		if(length<glueRTruncate) length = glueRTruncate;
		double[] dir = {(x2-x1)/length, (y2-y1)/length};
		double f = glueConstant*(length-glueR0);
		return new double[] {dir[0]*f,dir[1]*f};
	}

	// calculate LJ force between an edge and a node
	public void leonardJones() {
		// put every epithelial node in a tree structure
		blobTree = new PointQuadtree<BlobC>(blobCs.get(0),0,0,width,height);
		for(int i = 1; i < blobCs.size(); i++) {
			blobTree.insert(blobCs.get(i));
		}

		// for every wall segment, find the epithelial nodes that are close enough to exert LJ force
		for(int j = 0; j < wallCells.size(); j++) {
			Cell intraCell = wallCells.get(j);
			BlobC blob1, blob2;
			blob2 = intraCell.getBlob(0);
			for(int i = 0; i < intraCell.getSize(); i++) {
				blob1 = blob2;
				blob2 = intraCell.getBlob((i+1)%intraCell.getSize());
				List<BlobC> blobTrack = new ArrayList<>();
				// find the nodes
				blobTrack = blobTree.findInCircle((blob1.getX()+blob2.getX())/2, (blob1.getY()+blob2.getY())/2, rMax);
				// for every such node, apply LJ force
				for(BlobC blob : blobTrack) {
					if(blob == blob1 || blob == blob2) continue;
					double f[] = leonardJonesHelper(blob1,blob2,blob,leonardJonesConstant,leonardJonesR0);
					blob.applyLoad(f[0], f[1], drawLeonard);
				}
			}
		}

		// for every epithelial segment, find the epithelial nodes that are close enough to exert LJ force
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
	

	// calculate the LJ force between an edge and a node. The edge is defined by two nodes: blob1 & blob2. The node is blob 3. kij is LJ constant, r0 is resting length
	public double[] leonardJonesHelper(Blob blob1, Blob blob2, Blob blob3, double kij, double r0) {
		double x1 = blob1.getX(), y1 = blob1.getY();
		double x2 = blob2.getX(), y2 = blob2.getY();
		double x3 = blob3.getX(), y3 = blob3.getY();
		double rij = 0.0;
		// calculate the length of the edge
		double length = Math.sqrt(Math.pow(x2-x1,2)+Math.pow(y2-y1,2));
		// distance of the projection of node3 on the edge to node1 and node2
		double proj1 = (x3-x1)*(x2-x1)/length+(y3-y1)*(y2-y1)/length;
		double proj2 = (x3-x2)*(x2-x1)/length+(y3-y2)*(y2-y1)/length;
		double f = 0.0;
		double[] direction = {0.0, 0.0};
		if(Math.abs(proj1) > Math.abs(proj2)) {

			if (proj2 > 0) { // if the projection of node3 falls outside the edge
				rij = Math.sqrt((Math.pow(x3-x2, 2)+Math.pow(y3-y2, 2)));
				direction[0] = (x3-x2)/rij; direction[1] = (y3-y2)/rij;
			}else { // if the projection falls on the edge

				rij = Math.sqrt((Math.pow(x3-x1, 2)+Math.pow(y3-y1, 2)-Math.pow(proj1, 2)));
				if((x3-x1)*(y1-y2)+(y3-y1)*(x2-x1)>0){
					direction[0]=(y1-y2)/length; direction[1] = (x2-x1)/length;
				}else {
					direction[0]=(y2-y1)/length; direction[1] = (x1-x2)/length;
				}
			}
		}else {
			if(proj1 > 0) { // if the projection of node3 falls on the edge
				rij = Math.sqrt((Math.pow(x3-x1, 2)+Math.pow(y3-y1, 2)-Math.pow(proj1, 2)));
				if((x3-x1)*(y1-y2)+(y3-y1)*(x2-x1)>0){
					direction[0]=(y1-y2)/length; direction[1] = (x2-x1)/length;
				}else {
					direction[0]=(y2-y1)/length; direction[1] = (x1-x2)/length;
				}
			}else { // if the projection of node3 falls outside the edge
				rij = Math.sqrt((Math.pow(x3-x1, 2)+Math.pow(y3-y1, 2)));
				direction[0] = (x3-x1)/rij; direction[1] = (y3-y1)/rij;
			}
		}
		if(rij < rTruncate) rij = rTruncate;
//		f = kij/Math.pow(rij, 3)*(8*Math.pow(r0/rij, 6)-10*Math.pow(r0/rij, 3));
//		System.out.println("rij "+rij);

		f = kij/Math.pow(rij, 3)*(6*Math.pow(r0/rij, 4)-8*Math.pow(r0/rij, 2));
		double ff[] = {f*direction[0], f*direction[1]};
		return ff;
	}

	// not used
	// calculate the apical-basal elongation force on ectodermal cells
	public void ectoLoad() {
		// apical elongation ratio
		double i1 = ectoRatio[0];
		double i2 = ectoRatio[1];


		// basal elongation ratio
		double i3 = ectoRatio[2];
		double i4 = ectoRatio[3];
		// apical elongation spring constant
		double j1 = ectoConstant[0];
		double j2 = ectoConstant[1];
		double j1X = ectoConstant[0];
		double j2X = ectoConstant[1];

		//j1 = 0;
		//j2 = 0;
		j1 = j1/10;
		j2 = j2/10;



		// basal elongation spring constant
		double j3 = ectoConstant[2];
		double j4 = ectoConstant[3];
		System.out.println("ecto-"+j1+"-"+i1);

		double ectoRat;
		double ectoConst;


		int ectoNum = ectoEnd-ectoBegin;
		// calculate distribution of the ratio
		double func1[] = new double[ectoNum*nNode11];
		double func2[] = new double[ectoNum*nNode11];
		for(int i = 0; i < ectoNum*nNode11; i++) {
			func1[i] = (i2-i1)*(1.0/ectoNum/nNode11*i)+i1;
		}
		// calculate the distribution of the constant
		for(int j = 0; j < ectoNum*nNode11; j++) {
			func2[j] = (j2-j1)*(1.0/ectoNum/nNode11*j)+j1;
		}
		int count = 0;

		// apply load on apical side
		// apply load on the right half
		for(int i = 0; i < nMesh; i++) {
			Cell cell = cells.get(i);
			if(i <= ectoBegin)
			{
				System.out.println("Case 1\n");
				ectoRat = (i2-i1)*(1.0/ectoNum/nNode11*i)+i1;
				ectoConst = (j2X-j1X)*(1.0/ectoNum/nNode11*i)+j1X;
				for(int j = 0; j < nNode11; j++) {
					cell.ectoLoad(j,ectoConst,ectoRat);
					count++;
				}
			}
			else if (i >= nMesh-ectoBegin)
			{
				System.out.print("Case 2\n");
				ectoRat = (i2-i1)*(1.0/ectoNum/nNode11*i)+i1;
				ectoConst = (j2X-j1X)*(1.0/ectoNum/nNode11*i)+j1X;
				for(int j = 0; j < nNode11; j++) {
					cell.ectoLoad(j,ectoConst,ectoRat);
					count++;
				}
			}
			else
			{
				System.out.println("Case 3\n");
				ectoRat = (i2-i1)*(1.0/ectoNum/nNode11*i)+i1;
				ectoConst = (j2-j1)*(1.0/ectoNum/nNode11*i)+j1;
				for(int j = 0; j < nNode11; j++) {
					cell.ectoLoad(j,ectoConst,ectoRat);
					count++;
				}
			}
		}


		// apply load on basal side
		func1 = new double[ectoNum*nNode12];
		func2 = new double[ectoNum*nNode12];
		// apply load on the right half
		for(int i = 0; i < ectoNum*nNode12; i++) {
			func1[i] = (i4-i3)*(1.0/ectoNum/nNode12*i)+i3;
		}
		for(int j = 0; j < ectoNum*nNode12; j++) {
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
		// apply load on the left half
		for(int i = nMesh - ectoBegin-1; i >= nMesh - ectoEnd; i--) {
			Cell cell = cells.get(i);
			for(int j = nNode11+nNode2; j < nNode11 + nNode2+ nNode12; j++) {
				cell.ectoLoad(j,func2[count],func1[count]);
				count++;
			}
		}
	}
	
	// apply lateral load
	public void lateralLoad() {
		// lateral constriction ratio
		double i1 = lateralPara[0];
		double i2 = lateralPara[1];
//		double i1 = lateralRatio[0];
//		double i2 = lateralRatio[1];
		// lateral constriction spring constant
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
		for(int i = lateralBegin; i <lateralEnd; i++) {
			Cell cell = cells.get(i);
			for(int j = nNode11; j < nNode11+nNode2; j++) {
				cell.lateralLoad(j, func2[j-nNode11], func1[j-nNode11]);
			}
//			for(int j = nNode11+nNode12+nNode2; j < nNode11+2*nNode2+nNode12; j++) {
//				cell.compressiveLoad(j,func2[j-nNode11-nNode12-nNode2], func1[j-nNode11-nNode12-nNode2]);
//			}
		}
	}

	// apply apical constriction
	public void constrictLoad() {
		double i1 = constrictRatio[0];
		double i2 = constrictRatio[1];
		double j1 = constrictConstant[0];
		double j2 = constrictConstant[1];
		// time-dependent constriction constant
		j1 = constrictConstPara;
		System.out.println("compressPara-"+j1);
		double[] func1 = new double[constrictNum*nNode11];
		double[] func2 = new double[constrictNum*nNode11];
		for(int i = 0; i < constrictNum*nNode11; i++) {
			func1[i] = (i2-i1)*(1.0/constrictNum/nNode11*i)+i1;
		}
		double cutoff = j1* Math.exp(-0.5*j2*j2*(constrictNum*nNode11-1)/nNode11/constrictNum);
		for(int j = 0; j < constrictNum*nNode11; j++) {
//			func2[j] = j1*compressDistribution[j];
			double loc = j2/constrictNum/nNode11*j;
			
			func2[j] = j1* Math.exp(-0.5*loc*loc)-cutoff;
//			func2[j] = 0.75*func2[j];
		}
		

//		func2 is compress constant
//		func1 is compress ratio
		int count = 0;
		for(int i = 0; i<constrictNum; i++) {
			Cell cell = cells.get(i);
			for(int j = 0; j<nNode11; j++) {
				cell.constrictLoad(j, func2[count], func1[count]);
//				cell.compressiveLoad(j, compressConstant, compressRatio[i]);
				count++;
			}
		}
		count = 0;
		for(int i = cells.size()-1; i>=cells.size()-constrictNum; i--) {
			Cell cell = cells.get(i);
			for(int j = nNode11-1; j>=0; j--) {
				cell.constrictLoad(j, func2[count], func1[count]);
//				cell.compressiveLoad(j, compressConstant, compressRatio[cells.size()-1-i]);
				count++;
			}
		}
		
	}
	//Writer
	/**
	 * DrawingGUI method, here having all the cells take a step
	 */
	// write the location of the nodes to a file
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
			//long timestamp = System.nanoTime();
			String timestampString = String.valueOf(timeCount);
			WriteFileAtTimePoint.WriteToFile(writer, timestampString, s);
			/*try {
				writer.write(s);
				writer.newLine();
			}catch(IOException e1){
				e1.printStackTrace();
			}*/


		}
		if(timeCount==stopWriting) {
			try {
				writer.close();
			}catch(IOException e1) {
				e1.printStackTrace();
			}

		}
	}

	// this function is to control what happens when the timer ticks.
	// each time timer ticks, we call this function
	@Override
	public void handleTimer() {
		// calculate the value of time-dependent variables
		timeDependent();
		// apply all the mechanisms
		//internal cell forces:
		int count = 0;
		for(int i = 0; i<constrictNum; i++) {
			Cell cell = cells.get(i);
			for(int j = 0; j<nNode11; j++) {
				cell.elastic(elasticConstant/1);
				//cell.cornerstiff(stiffnessConstant/0.7, nNode11, nNode12, nNode2);
				cell.hasLoad = true;
//				cell.compressiveLoad(j, compressConstant, compressRatio[i]);
				count++;
			}
		}
		count = 0;
		for(int i = cells.size()-1; i>=cells.size()-constrictNum; i--) {
			Cell cell = cells.get(i);
			for(int j = nNode11-1; j>=0; j--) {
				cell.elastic(elasticConstant/1);
				//cell.cornerstiff(stiffnessConstant/0.7, nNode11, nNode12, nNode2);
				cell.hasLoad = true;
//				cell.compressiveLoad(j, compressConstant, compressRatio[cells.size()-1-i]);
				count++;
			}
		}
		for (Cell cell : cells) {
			//elastic force
			//cell.elastic(elasticConstant);
			if(!cell.hasLoad) {
				cell.elastic(elasticConstant);
				//cell.cornerstiff(stiffnessConstant, nNode11, nNode12, nNode2);
			}else
			{cell.hasLoad = false;}

			cell.osmosis(osmosisConstant);	//volume conservation
			cell.cornerstiff(stiffnessConstant, nNode11, nNode12, nNode2);
			
		}

		//external cell forces (cell-cell):
		leonardJones();

		//cause invagination
		constrictLoad();
		lateralLoad();

//		apicalGlue();
		//spring constant applied here
		ectoLoad();		//from ectodermal cells

//		yolkCell.osmosis(yolkOsmosisConstant);
		yolkCell.hydrostatic(yolkHydroConstant);	//hydrostatic (in yolk internal areas enclosed by cell)

		// ask the nodes to move
		for (BlobC blob : blobCs) {
			blob.step(delay/10000.0);
		}
		writeFile();
		// Now update the GUI.
		repaint();

	}

	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				new CellsGUI();
			}
		});
	}
}

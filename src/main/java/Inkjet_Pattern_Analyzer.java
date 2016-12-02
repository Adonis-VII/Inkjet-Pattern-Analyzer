/*
Inkjet_Pattern_Analyzer is an ImageJ plugin that analyzes optical images of calibration grids
and determines the reliability of the printing process.
Copyright (C) 2015  Ryan P. Sullivan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.*;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.BackgroundSubtracter;
import ij.process.*;
import java.awt.Color;
import java.awt.Font;
import java.awt.Label;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import ij.process.AutoThresholder.Method;
import ij.process.ImageConverter;
import ij.util.Tools;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;

public class Inkjet_Pattern_Analyzer implements PlugIn {

int rows, columns;
int countExpected, countActual, countSortedData; //Tracks expected number of features (rows * columns), the found number of features, and sorted number of features + buffers
double radius;
double distance, ActiveDistance;
double deviations = 3;
double[][] Data, SortedData, OrderData;
protected ImagePlus Image, ActiveImage, BinaryImage;
ImageProcessor ActiveImageProcessor;
ImageConverter ActiveImageConverter;
ResultsTable rt, rtPassthrough;
double minSize, minSize_Default = 200.0, maxSize, maxSize_Default = Double.POSITIVE_INFINITY;
boolean intelligent, binary, outline, showResults, batch;	//Variables for dialog window
boolean logging = true, Progresslogging = true;	
String[] BackgroundRemovalMethods = {"None","Manual","Automatic","Manual+Automatic"};
String[] ThresholdingMethods = {"None","Otsu"};
String[] SortingMethods = {"Insertion","Grid","Neighbor"};
String BackgroundRemovalMethod, ThresholdingMethod, SortingMethod;

String InputDirectory, OutputDirectory;
String[] imagelist;
String imageName;
int filenum;

	//First method called by ImageJ when run.  Main() isn't used in ImageJ plugins.
	public void run(String arg) {
		if(!showSetup()){
			return;
		}
		if(batch) {
			setDirectories();
			for(int i=0; i<filenum; i++) {
				getImage(i);
				if(!ProcessImage()){
					IJ.log("Error: Could not analyze..."+imageName);	
				}
				else{
					//CalculateOrder();
					writeOutputs(i);
					IJ.log(imageName + " Processed");
				}
			}
		}
		else {
			Image=IJ.getImage();
			imageName = Image.getTitle();
			reset();
			if(imageName.endsWith(".tiff")){
			imageName = imageName.substring(0,imageName.indexOf(".tiff"));
			}
			if(imageName.endsWith(".tif")){
				imageName = imageName.substring(0,imageName.indexOf(".tif"));
			}
			if(Progresslogging){
				IJ.log("Working on " +imageName);
			}
			if(!ProcessImage()) {
				IJ.log("Error: Could not analyze");	
				return;
			}
			displayOutputs();	
		}
		IJ.log("Finished");
		return;	
	}
	//Setup method that displays a GUI and asks for parameters
	public boolean showSetup() {
		GenericDialog gd = new GenericDialog("Grid Analyzer (c) 2015 R. Sullivan");
		gd.addChoice("Background Removal: ",BackgroundRemovalMethods,BackgroundRemovalMethods[0]);
		gd.addChoice("Thresholding: ",ThresholdingMethods,ThresholdingMethods[1]);
		gd.addChoice("Sorthing Method: ",SortingMethods,SortingMethods[0]);
		gd.addStringField("Size (pixel^2):", minSize_Default+"-"+maxSize_Default,10);
		gd.addNumericField("Colums:",6,0);
		gd.addNumericField("Rows:",5,0);
		gd.addNumericField("Diameter",15,0);
		gd.addNumericField("Distance",80,0);
		gd.addCheckbox("Intelligent?", true);
		gd.addCheckbox("Show Binary?", false);
		gd.addCheckbox("Show Outlines?", true);
		gd.addCheckbox("Show Results?", true);
		gd.addCheckbox("Batch Process?", false);
		gd.showDialog();
		if (gd.wasCanceled()){
			return false;
		}
		if (gd.invalidNumber()) {
			IJ.error("Invalid Number!");
			return false;
		}
		BackgroundRemovalMethod = gd.getNextChoice();
		ThresholdingMethod = gd.getNextChoice();
		SortingMethod = gd.getNextChoice();
		columns = (int)gd.getNextNumber();
		rows = (int)gd.getNextNumber();
		Double DiameterInput= gd.getNextNumber();
		distance = gd.getNextNumber();
		ActiveDistance = distance;
		//min-max
		String size = gd.getNextString();
		String[] MinMax = Tools.split(size,"-");
		double min = (MinMax.length >= 1) ? gd.parseDouble(MinMax[0]) : 0.0;
		double max = (MinMax.length == 2) ? gd.parseDouble(MinMax[1]) : Double.NaN;
		minSize = min;
		maxSize = max;
		//minSize = Double.isNaN(min) ? minSize_Default : mins/unitSquared;
        	//maxSize = Double.isNaN(max) ? maxSize_Default : maxs/unitSquared;
        	//if (minSize<DEFAULT_MIN_SIZE) minSize = DEFAULT_MIN_SIZE;
        	//if (maxSize<minSize) maxSize = maxSize_Default;
        radius = Math.round(DiameterInput.intValue()/2*1.5);
		countExpected = rows*columns;
		intelligent = gd.getNextBoolean();
		binary = gd.getNextBoolean();
		outline = gd.getNextBoolean();
		showResults = gd.getNextBoolean();
		batch = gd.getNextBoolean();
		return true;
	}
	//Prompts user to choose input and output directories, and creates a list of all valid .tiff files in the input directory
	public void setDirectories() {
		DirectoryChooser directory = new DirectoryChooser("Choose input directory");
		InputDirectory = directory.getDirectory();
		directory = new DirectoryChooser("Choose output directory");
		OutputDirectory = directory.getDirectory();
		File folder = new File(InputDirectory);
		File[] files = folder.listFiles(new FilenameFilter() {
			@Override
            	public boolean accept(File dir, String name) {
            		//ignore outline.tiff & outline.tif
                	if((name.toLowerCase().endsWith(".tiff") && !name.toLowerCase().endsWith("outline.tiff")) || (name.toLowerCase().endsWith(".tif") && !name.toLowerCase().endsWith("outline.tif"))){
                		return true;
                	} else {
                    	return false;
                	}
            	}
        	});
        	filenum=files.length;
        	imagelist = new String[filenum];
        	for(int i=0; i<filenum; i++)
        		imagelist[i]=files[i].getName();
		return;
	}
	//Reset ActiveImage to current Image
	public void reset() { 
		ActiveImage=Image.duplicate();
		ActiveImageProcessor = ActiveImage.getProcessor();
		ActiveImageConverter = new ImageConverter(ActiveImage);
		return;
	}
	//Sets Image + global parameters based on directory
	public void getImage(int index) {
		Opener opener = new Opener();
		Image = opener.openImage(InputDirectory + imagelist[index]);
		reset();
		imageName = Image.getTitle();
		if(imageName.endsWith(".tiff")){
			imageName = imageName.substring(0,imageName.indexOf(".tiff"));
		}
		if(imageName.endsWith(".tif")){
			imageName = imageName.substring(0,imageName.indexOf(".tif"));
		}
		if(Progresslogging){
			IJ.log("Working on " + imageName);
		}
		return;
	}
	//Main method to process image
	public boolean ProcessImage(){
		if (BackgroundRemovalMethod=="Manual" || BackgroundRemovalMethod=="Manual+Automatic"){
			BackgroundSubtracter subtracter = new BackgroundSubtracter();
			subtracter.rollingBallBackground(ActiveImageProcessor, radius, false, true,true,false,true);
			if(Progresslogging)
				IJ.log("Removed Background with radius " + radius);
		}
		if (BackgroundRemovalMethod=="Automatic"|| BackgroundRemovalMethod=="Manual+Automatic") {//interative //check if area is * 10 larger than average
            ActiveImageConverter.convertToGray16();
            IJ.run(ActiveImage, "Auto Threshold", "method=Otsu white");
            //ActiveImageProcessor.setAutoThreshold(AutoThresholder.Method.valueOf("Otsu"), true);
            if(!setMeasurementArray()){
            	return false;
            }
			rtPassthrough = rt;
			rtPassthrough.show("Results2");
			int count = StDevFilter(6);
			IJ.log("Removed " + count + " Feret Outliers, " + countActual + " features left.");
			radius = CalculateMax(Data[6])*1.5;
			reset();
			BackgroundSubtracter subtracter = new BackgroundSubtracter();
			subtracter.rollingBallBackground(ActiveImageProcessor, radius, false, true,true,false,true);
			if(Progresslogging)
				IJ.log("Removed Background with radius " + radius);
		}
		if(ThresholdingMethod=="Otsu"){
            ActiveImageConverter.convertToGray16();
            IJ.run(ActiveImage, "Auto Threshold", "method=Otsu white");
            //ActiveImage.getProcessor().setAutoThreshold(AutoThresholder.Method.valueOf("Otsu"), true);
            //ActiveImage.updateImage();
            BinaryImage = ActiveImage.duplicate();
            IJ.log("Thresholded");
		}
		if(!setMeasurementArray()){
            return false;
        }
		//check number of expected
		if(countActual>countExpected) { //if excess number of dots try to filter by area
			StDevFilter(0);
			if(countActual>countExpected){
				IJ.log("Error: Unexpected excess of dots.");	
				return false;
			}
		}
		if(!intelligent) {//Do exactly as user specifices
			if(SortingMethod=="Insertion"){
				InsertionSort();
			}
			else if(SortingMethod=="Grid"){
				GridSort();
			}
			else{
				NeighborSort();
			}
		}
		else{//Try and figure optimal sorting
			if(countActual==countExpected) {
				if(SortingMethod=="Insertion"){
					InsertionSort();
				}
				else if(SortingMethod=="Grid"){
					GridSort();
				}
			}
			//consider checking for dots twice the volume
			else if(countActual < countExpected) {
				GridSort();  //used for setting distances
				removeDistanceOutliers();
				if (!GridSort())
					return false;
			}
		}
		CalculateOrder();
		return true;
	}
	//Display any output windows
	public void displayOutputs() {
		if(outline){
			Outline();
			ActiveImage.setTitle(imageName + "_outline");
			ActiveImage.show();
		}
		if(binary){
			BinaryImage.setTitle(imageName + "_binary");
			BinaryImage.show();
		}
		if(showResults){
			OutputToResults();
			rt.show("Results");
		}
		//rtPassthrough.show("Results Passthrough");
		return;
	}
	//Writes any output windows to disk
	public void writeOutputs(int index) {
		if(outline){
			Outline();
			FileSaver saver = new FileSaver(ActiveImage);
			saver.saveAsTiff(OutputDirectory+imageName+"_Outline.tiff");
		}
		if(showResults){
			OutputToResults();
			try{
				rt.saveAs(OutputDirectory+imageName+"_Results.csv");
			}catch(IOException e){
				IJ.log("Could not write results for "+imageName);
			}
		}
		return;
	}
	//Insertion Sort of grid with expected number of dots
	public void InsertionSort() {
		countSortedData = countActual;
		SortedData = Data;
		for(int i=0; i<rows; i++) {//step by row
			for(int j=i*columns+1; j<columns*(i+1); j++) {//insertion sort
				double areatemp = SortedData[0][j];
				double xtemp=SortedData[1][j];
				double ytemp=SortedData[2][j];
				double perimetertemp = SortedData[3][j];
				double circularitytemp = SortedData[4][j];
				double roundnesstemp = SortedData[5][j];
				double FeretTemp = SortedData[6][j];
				int k;
				for(k=j-1; (k>=columns*i) && (SortedData[1][k] > xtemp); k--) {
					SortedData[0][k+1] = SortedData[0][k];
					SortedData[1][k+1] = SortedData[1][k];
					SortedData[2][k+1] = SortedData[2][k];
					SortedData[3][k+1] = SortedData[3][k];
					SortedData[4][k+1] = SortedData[4][k];
					SortedData[5][k+1] = SortedData[5][k];
					SortedData[6][k+1] = SortedData[6][k];
				}
				SortedData[0][k+1]=areatemp;
				SortedData[1][k+1]=xtemp;
				SortedData[2][k+1]=ytemp;
				SortedData[3][k+1]=perimetertemp;
				SortedData[4][k+1]=circularitytemp;
				SortedData[5][k+1]=roundnesstemp;
				SortedData[6][k+1]=FeretTemp;
			}			
		}
		if(Progresslogging)
			IJ.log("Sorted grid using Insertion Sort");
		return;
	}
	//Neighbor Sort of grid
	public boolean NeighborSort() {
		if(Progresslogging)
			IJ.log("Attempting Neighbor Sort...");
		int[][] RelativeCoordinates = new int[2][countActual];
		RelativeCoordinates [0][0]=0; //Set initial dot x pos at 0
		RelativeCoordinates [1][0]=0; //Set initial dot y pos at 0
		int[][] relativeNeighborPosition = {{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{-1,1},{0,1},{1,1}};
		int errorScalar = 2;
		for(int i =0; i<countActual; i++) {//step through each dot
			double[][] idealNeighbors = new double[2][8];
			List<Integer> neighborList = new ArrayList<Integer>();
			for(int j =0; j<8; j++){//set ideal neighbor positions
				idealNeighbors[0][j] = Data[1][j] + ActiveDistance*relativeNeighborPosition[j][0];
				idealNeighbors[1][j] = Data[2][j] + ActiveDistance*relativeNeighborPosition[j][1];
			}
			for(int j = 0; j<i; j++) {//find all neighbors within distance*errorScalar (2 by default) 
				if(Data[1][j]>Data[1][i]-ActiveDistance*errorScalar &&  Data[1][j]<Data[1][i]-ActiveDistance*errorScalar){
					if(Data[2][j]>Data[2][i]-ActiveDistance*errorScalar &&  Data[2][j]<Data[2][i]-ActiveDistance*errorScalar){
						neighborList.add(j);
					}
				}
			}
			int index = 0;
			for(int j=0; j<neighborList.size(); j++) {//compare neighbors actual position to ideal neighbors 0,6,7
				int num = neighborList.get(j);
				double difference = CalculateDistance(Data[1][num],Data[2][num],idealNeighbors[0][0],idealNeighbors[1][0]);
				for(int k=1; k<8; k++) {//compare given neighbor to ideal positions, finding the closest match
					double euclidianDistance=CalculateDistance(Data[1][num],Data[2][num],idealNeighbors[0][k],idealNeighbors[1][k]);
					if(euclidianDistance<difference) {
						difference = euclidianDistance;
						index = k;
					}
				}
				//set relative coordinate of new dot
				RelativeCoordinates[0][j]=RelativeCoordinates[0][index]+relativeNeighborPosition[index][0];
				RelativeCoordinates[1][j]=RelativeCoordinates[1][index]+relativeNeighborPosition[index][1];
			}
		}
		return false;
	}
	public boolean GridSort() { //used for initial alignment to use when calculating distances for second run.
		if(Progresslogging)
			IJ.log("Attempting Grid Sort...");
		int[][] RelativeCoordinates = new int[2][countActual];
		RelativeCoordinates [0][0]=0; //Set initial dot x pos at 0
		RelativeCoordinates [1][0]=0; //Set initial dot y pos at 0
		int minXShift = 0; //min relative X position
		int maxXShift = 0; //max relative X position
		int minYShift = 0; //min relative Y position
		int maxYShift = 0; //max relative Y position
		for(int i =1; i<countActual; i++) { //constant time O(1)
			double xReal = Data[1][i]-Data[1][0];  //Calculates x-distance given dot and initial dot
			double yReal = Data[2][i]-Data[2][0];  //Calculates y-distance given dot and initial dot
			int xGrid = (int) Math.round(xReal/ActiveDistance);  //Divides by expected distance and rounds to integer
			int yGrid = (int) Math.round(yReal/ActiveDistance);  //Divides by expected distance and rounds to integer
			RelativeCoordinates[0][i]=xGrid;	//Sets x-grid offset for given drop to array
			RelativeCoordinates[1][i]=yGrid;	//Sets y-grid offset for given drop to array
			if (minXShift > xGrid) minXShift = xGrid;
			if (maxXShift < xGrid) maxXShift = xGrid;
			if (minYShift > yGrid) minYShift = yGrid;
			if (maxYShift < yGrid) maxYShift = yGrid;
		}
		int maxXCount = maxXShift-minXShift+1;
		int maxYCount = maxYShift-minYShift+1;
		IJ.log("Grid of " + maxXCount + " x " + maxYCount + " formed");
		if(minXShift < 0 || minYShift <0) {//shift grid so top row, and left column are 0th
			IJ.log("Shift!");
			for(int i =0; i<countActual;i++) {
				if (minXShift < 0){
					RelativeCoordinates[0][i] = RelativeCoordinates[0][i] - minXShift;
				}
				if (minYShift < 0){
					RelativeCoordinates[1][i] = RelativeCoordinates[1][i] - minYShift; 
				}
			}
		}
		countSortedData = (maxXCount*maxYCount > countActual) ? maxXCount*maxYCount:countActual;
		SortedData = new double[7][countSortedData];
		int count=0;
		for(int i =0; i<countActual; i++) {
			int position = columns*RelativeCoordinates[1][i]+RelativeCoordinates[0][i];
			SortedData[0][position]=Data[0][i]; //Set Area
			if(SortedData[0][position] !=0){
				if(count<position){
					count=position;
				}
			}
			SortedData[1][position]=Data[1][i]; //Set X
			SortedData[2][position]=Data[2][i]; //Set Y
			SortedData[3][position]=Data[3][i]; //Set Perimeter
			SortedData[4][position]=Data[4][i]; //Set Circularity
			SortedData[5][position]=Data[5][i]; //Set Roundness
			SortedData[6][position]=Data[6][i]; //Set Feret
		}
		count++;
		//countSortedData = count;
		if(Progresslogging){
			IJ.log(countActual + " features mapped to " + countSortedData +" places.");
		}
		return true;
	}
	//finds distances Depricated
	public void findDistance() {
		int length=(countActual*(countActual-1))/2;
		int[] Lengths = new int[length];
		int count =0;
		for(int i =0; i < countActual-1; i++) {
			for(int j = i+1; j<countActual; j++) {
				int dist =(int) Math.pow(Math.pow(Data[2][j]-Data[2][i],2) + Math.pow(Data[1][j] -Data[1][i],2),0.5);
				Lengths[count]=dist;
				count++;
			}
		}		
		Arrays.sort(Lengths);
		for(int i =0;i<count;i++)
			IJ.log(""+Lengths[i]);
	}
	//filters outliers from an array
	public int StDevFilter(int index) {//need to add interation in cases where count isn't 0
		if(Progresslogging){
			if(index == 0)
				IJ.log("Attempting to remove Area Outliers...");
			else if(index == 6)
				IJ.log("Attempting to remove Feret Outliers...");
		}
		double[] tempArray = Data[index];
		int count = 0;
		double mean = CalculateMean(tempArray);
		double stDev = CalculateStandardDev(tempArray, mean);
		double[] tempArrayFilter = new double[countActual];  
		for(int i=0; i<countActual; i++) { //sets up temperorary array 1 or 0 values to act as a mask for the next loop
			if(Math.abs(tempArray[i]-mean) > stDev*deviations){
				tempArrayFilter[i]=1;
				count++;
			}
		}
		if(count>0){
			double[][] tempData = new double[7][countActual-count];
			int j =0;
			for(int i=0; i<countActual; i++) {
				if(tempArrayFilter[i] == 0) {
					tempData[0][j]=	Data[0][i];
					tempData[1][j]=	Data[1][i];
					tempData[2][j]=	Data[2][i];
					tempData[3][j]=	Data[3][i];
					tempData[4][j]=	Data[4][i];
					tempData[5][j]=	Data[5][i];
					tempData[6][j]= Data[6][i];
					j++;
				}
			}
			Data = new double[7][countActual-count];
			Data = tempData;
		}
		countActual = countActual-count;
		return count;
	}
	//calculates the median value of an array
	public int CalculateMedian(double[] data){
		double Median;
		Arrays.sort(data);
		int dataCount = data.length;
		if(dataCount % 2 ==0) {//even
			int num = (int) (dataCount/2 + .5);
			Median = (data[num]+data[num+1])/2;
		}
		else {
			Median = data[dataCount/2];
		}
		int median = (int) Median;
		return median;
	}
	//uses sortedData to determine distance distribution
	public void removeDistanceOutliers() {
		if(Progresslogging)
			IJ.log("Attempting to remove Distance Outliers...");
		double[][] DistanceData = new double[2][countSortedData];
		DistanceData[0]=CalculateX();
		DistanceData[1]=CalculateY();
		double[] DistanceData2 = new double[countSortedData*2];
		System.arraycopy(DistanceData[0],0,DistanceData2,0,countSortedData-1);
		System.arraycopy(DistanceData[1],0,DistanceData2,countSortedData,countSortedData-1);
		ActiveDistance = CalculateMean(DistanceData2);
		//DecimalFormat df = new DecimalFormat("#.00");
		if(Progresslogging)
			IJ.log("Distance: " + distance + " changed to " + ActiveDistance);
		return;
	}
	//Grabs data from ParticleAnalyzer's result table and sets it as an array
	public boolean setMeasurementArray() {
		rt = new ResultsTable();
		int options = ParticleAnalyzer.CLEAR_WORKSHEET;
		int measurements = Measurements.AREA + Measurements.CENTROID + Measurements.PERIMETER +  Measurements.SHAPE_DESCRIPTORS + Measurements.FERET;
		ParticleAnalyzer ActiveImageAnalyzer = new ParticleAnalyzer(options,measurements,rt,minSize,maxSize);
		if(!ActiveImageAnalyzer.analyze(ActiveImage)){
			IJ.log("Could not grab Measurements");
			return false;
		}
		countActual = rt.size();
		Data = new double[7][countActual];
		Data[0] = rt.getColumnAsDoubles(rt.getColumnIndex("Area"));//set Area
		Data[1] = rt.getColumnAsDoubles(rt.getColumnIndex("X"));//set X
		Data[2] = rt.getColumnAsDoubles(rt.getColumnIndex("Y"));//set Y
		Data[3] = rt.getColumnAsDoubles(rt.getColumnIndex("Perim.")); //set Perimeter
		Data[4] = rt.getColumnAsDoubles(rt.getColumnIndex("Circ."));//set Circularity
		Data[5] = rt.getColumnAsDoubles(rt.getColumnIndex("Round"));//set Roundness
		Data[6] = rt.getColumnAsDoubles(rt.getColumnIndex("Feret"));//set Feret
		IJ.log(countActual + " features found.");
		return true;
	}
	//Outlines and labels current image
	public void Outline() {
		if(Progresslogging)
			IJ.log("Outlining...");
		IJ.run(ActiveImage, "Outline", "");
		ActiveImageProcessor = ActiveImage.getProcessor();
		for(int i=0; i<countSortedData;i++) {
			if(SortedData[0][i] != 0){
				String stamp = Integer.toString(i+1);
				int x = (int)Math.round(SortedData[1][i] - ActiveImageProcessor.getStringWidth(stamp)/2);
				Font font = ActiveImageProcessor.getFont();
				int y = (int)Math.round(SortedData[2][i]+ActiveImageProcessor.getFontMetrics().getHeight()/2);
				ActiveImageProcessor.drawString(stamp, x,y,Color.white);
			}
		}
	}
	//Writes the final data to Result Table using ImageJ's methods
	public void OutputToResults() {
		if(Progresslogging)
			IJ.log("Outputing Results..");
		rt = new ResultsTable();
		for(int i = 0; i<countSortedData; i++) {
			rt.incrementCounter();
			rt.addValue("Area", SortedData[0][i]);
			if(OrderData[0][i] != 0)
				rt.addValue("Theta", OrderData[0][i]);
			else
				rt.addValue("Theta", "");
			if(OrderData[1][i] !=0)
				rt.addValue("Delta X",OrderData[1][i]);
			else
				rt.addValue("Delta X", "");
			if(OrderData[2][i] !=0)
				rt.addValue("Delta Y",OrderData[2][i]);
			else
				rt.addValue("Delta Y", "");
			rt.addValue("Perimeter",SortedData[3][i]);
			rt.addValue("Circularity",SortedData[4][i]);
			rt.addValue("Roundness",SortedData[5][i]);
			rt.addValue("     "," ");
			rt.addValue(" ", " ");
			rt.addValue("Average", " ");
			rt.addValue("St Dev", " ");
		}
		double AreaAverage = CalculateMean(SortedData[0]);
		double AreaStDev = CalculateStandardDev(SortedData[0],AreaAverage);
		rt.setValue(8, 0, "Area");
		rt.setValue(9, 0, AreaAverage);
		rt.setValue(10, 0, AreaStDev);
		double ThetaAverage = CalculateMean(OrderData[0]);
		double ThetaStDev = CalculateStandardDev(OrderData[0],ThetaAverage);
		rt.setValue(8, 1, "Theta");
		rt.setValue(9, 1, ThetaAverage);
		rt.setValue(10, 1, ThetaStDev);
		double XAverage = CalculateMean(OrderData[1]);
		double XStDev = CalculateStandardDev(OrderData[1],XAverage);
		rt.setValue(8, 2, "X");
		rt.setValue(9, 2, XAverage);
		rt.setValue(10, 2, XStDev);
		double YAverage = CalculateMean(OrderData[2]);
		double YStDev = CalculateStandardDev(OrderData[2],YAverage);
		rt.setValue(8, 3, "Y");
		rt.setValue(9, 3, YAverage);
		rt.setValue(10, 3, YStDev);	
		double PerimeterAverage = CalculateMean(SortedData[3]);
		double PerimeterStDev = CalculateStandardDev(SortedData[3],PerimeterAverage);
		rt.setValue(8, 4, "Perimeter");
		rt.setValue(9, 4, PerimeterAverage);
		rt.setValue(10, 4, PerimeterStDev);	
		double CircularityAverage = CalculateMean(SortedData[4]);
		double CircularityStDev = CalculateStandardDev(SortedData[4],CircularityAverage);
		rt.setValue(8, 5, "Circularity");
		rt.setValue(9, 5, CircularityAverage);
		rt.setValue(10, 5, CircularityStDev);	
		double RoundnessAverage = CalculateMean(SortedData[5]);
		double RoundnessStDev = CalculateStandardDev(SortedData[5],RoundnessAverage);
		rt.setValue(8, 6, "Roundness");
		rt.setValue(9, 6, RoundnessAverage);
		rt.setValue(10, 6, RoundnessStDev);	
		//Parameters
		rt.setValue(8,8,"Background Removal");
		rt.setValue(9,8,BackgroundRemovalMethod);
		rt.setValue(8,9, "Radius");
		if(BackgroundRemovalMethod=="None"){
			rt.setValue(9,9,"N/A");
		}
		else{
			rt.setValue(9,9,""+radius);
		}
		rt.setValue(8,10,"Sorting Method");
		rt.setValue(9,10,"Grid Sort");
	}
	//Calculates the geometric order values and saves them to OrderData array
	public void CalculateOrder() {
		if(Progresslogging)
			IJ.log("Calculating Order Values...");
		OrderData = new double[3][countSortedData];
		OrderData[0]=CalculateTheta();
		OrderData[1]=CalculateX();
		OrderData[2]=CalculateY();
	}
	//Calculates the theta values for a given sorted array of features
	public double[] CalculateTheta() {
		double theta[] = new double[countSortedData];
		for(int i =0; i<rows-1;i++) { //set row
			for(int j=i*columns;j<(i+1)*columns-1;j++) { //set dot
				if(SortedData[0][j] != 0 && SortedData[0][j+1] !=0 && SortedData[0][j+columns] != 0) {
					double angle1 = Math.atan((SortedData[2][j+1]-SortedData[2][j])/(SortedData[1][j+1]-SortedData[1][j]));
					double angle2 = Math.atan((SortedData[2][j]-SortedData[2][j+columns])/(SortedData[1][j]-SortedData[1][j+columns]));
					double result = Math.abs((angle1+angle2)*57.2958);
					//double check
					theta[j]= result;
				}
			}
		}
		return theta;
	}
	//Calculates the deltaX values for a given sorted array of features
	public double[] CalculateX() {
		double x[] = new double[countSortedData];
		for(int i =0; i<rows;i++) { //set row
			for(int j=i*columns; j<(i+1)*columns-1;j++) {//set dot
				if(SortedData[0][j] != 0 && SortedData[0][j+1] != 0) 
					x[j] = SortedData[1][j+1] - SortedData[1][j];
			}
		}
		return x;
	}
	//Calculates the deltaY values for a given sorted array of features
	public double[] CalculateY() {
		double y[] = new double[countSortedData];
		for(int i=0;i<rows-1;i++) {//set row
			for(int j=0; j<(i+1)*columns;j++) {//set dot
				if(SortedData[0][j] != 0 && SortedData[0][j+columns] !=0)
					y[j] = SortedData[2][j+columns] - SortedData[2][j];
			}
		}
		return y;

	}
	//Calculates the Euclidian distance between two points given their coordinates
	public double CalculateDistance(double x1, double y1, double x2, double y2) {
		double euclidianDistance = Math.sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
		return euclidianDistance;
	}
	//Calculates the Euclidian distance between two points given their indices
	public double CalculateDistance(int i, int j){
		double euclidianDistance = Math.sqrt((Data[1][j]-Data[1][i])*(Data[1][j]-Data[1][i])+(Data[2][j]-Data[2][i])*(Data[2][j]-Data[2][i]));
		return euclidianDistance;
	}
	//Calculates the mean of an array
	public double CalculateMean(double[] data) { //ignores 0 values
		double sum = 0, mean = 0;
		int count = 0;
		for(int i = 0; i<data.length; i++){
			if(data[i] != 0){
				sum = sum + data[i];
				count++;	
			}
		}
		mean = sum/count;
		return mean;
	}
	//Calculates the standard deviation of an array given the mean
	public double CalculateStandardDev(double[] data, double mean) { //ignores 0 values
		double sum = 0, stDev = 0;
		int count = 0;
		for(int i = 0; i<data.length; i++) {
			if(data[i] != 0){	
				sum = sum + Math.pow(data[i]-mean,2);
				count++;
			}
		}
		stDev = Math.sqrt(sum/count);
		return stDev;
	}
	//Finds the maximum value of an array
	public double CalculateMax(double[] data) {
		double max =0;
		for(int i =0; i<countActual; i++) {
			if(data[i]>max)
				max = data[i];
		}
		return max;
	}
}
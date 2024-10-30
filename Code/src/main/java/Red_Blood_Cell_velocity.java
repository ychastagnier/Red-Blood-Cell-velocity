import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Menus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.MessageDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.io.DirectoryChooser;
import ij.io.Opener;
import ij.measure.ResultsTable;

import java.awt.Button;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.io.IOException;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.ToolTipManager;

import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.plugin.Zoom;
import ij.plugin.frame.RoiManager;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.plugin.Colors;
import ij.plugin.Concatenator;
import ij.plugin.ImageCalculator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.stream.IntStream;


/**
 * Red Blood Cell velocity v1.0.0
 * @author Yan Chastagnier
 */
public class Red_Blood_Cell_velocity implements PlugIn {
	public void run(String arg0) {
		ToolTipManager.sharedInstance().setDismissDelay(30000);
		RBCVProcess.getInstance().execSubFunction(arg0);
	}
	
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Red_Blood_Cell_velocity.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);
		
		// create the ImageJ application context with all available services
		//final ImageJ ij = new ImageJ();
		new ImageJ();
		
		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
		//ij.quit();
	}
}

class RBCVProcess implements ActionListener, ItemListener, TextListener, WindowListener, MouseListener, FocusListener {
	private static RBCVProcess instance = null;
	private Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
	private double screenWidth = screenSize.getWidth();
	private double screenHeight = screenSize.getHeight();
	private double framePosX = Prefs.get("RBCv.framePosX", 0);
	private double framePosY = Prefs.get("RBCv.framePosY", 0);
	private double verifPosX = Prefs.get("RBCv.verifPosX", 0);
	private double verifPosY = Prefs.get("RBCv.verifPosY", 0);
	private Frame frame = null;
	private Frame verif = null;
	private boolean frameInitialized = false;
	private boolean verifInitialized = false;
	private Button selectDirectoryBtn = new Button("Select directory");
	private Choice depthValue = new Choice();
	private JLabel stackSuffixLabel = new JLabel("Stacks suffix ", JLabel.RIGHT);
	private TextField stackSuffixValue = new TextField(Prefs.get("RBCv.stackSuffix", ".tif"));
	private boolean stackSuffixChanged = false;
	//private JLabel shiftSuffixLabel = new JLabel(" Shifts suffix ", JLabel.RIGHT);
	//private TextField shiftSuffixValue = new TextField(Prefs.get("RBCv.shiftSuffix", "-shifts.txt"));
	//private boolean shiftSuffixChanged = false;
	private JLabel directoryLabel = new JLabel("Directory path ", JLabel.RIGHT);
	private TextField directoryValue = new TextField("directory path will be displayed here");
	private JLabel stacksListLabel = new JLabel("Stacks path ", JLabel.RIGHT);
	private Choice stacksListValue = new Choice();
	
	private Button openStackBtn = new Button("Open stack");
	private Button saveROIsBtn = new Button("Save ROIs");
	private Button exportROIsBtn = new Button("Export ROIs");
	private Button xtGenBtn = new Button("XT Gen");
	private Button measureVelocityBtn = new Button("Measure velocity");
	private Button measureVelocitySlidingBtn = new Button("Measure sliding velocity");
	private Button manualVerifBtn = new Button("Manual verif");
	private Button evaluateVerifBtn = new Button("Evaluate verif");
	private Button automaticVerifBtn = new Button("Auto verif");
	private Button rollingShutterBtn = new Button("Rolling shutter");
	private Button velocityLimitBtn = new Button("Velocity limit");
	private Button drawVelocityBtn = new Button("Draw velocity");
	
	private Button updateAngleBtn = new Button("Update angle [u]");
	private Button discardAngleBtn = new Button("Remove angle [r]");
	private Button openStackVerifBtn = new Button("Open stack");
	private Button openXTRawVerifBtn = new Button("Open raw XT");
	private Checkbox openVirtualCB = new Checkbox("Virtual stack", true);
	private Checkbox autoCloseCB = new Checkbox("Autoclose", true);
	private Button validateROIBtn = new Button("Validate ROI [y]");
	private Button unvalidateROIBtn = new Button("Discard ROI [n]");

	public final String VELOCITY_FILENAME = "ResultsVelocity.csv";
	public final String VELOCITY_MANUAL_FILENAME = "ResultsVelocityManual.csv";
	public final String VELOCITY_AUTO_FILENAME = "ResultsVelocityAuto.csv";
	public final String VELOCITY_SLIDING_FILENAME = "ResultsVelocitySliding.csv";
	public final String VELOCITY_SLIDING_AUTO_FILENAME = "ResultsVelocitySlidingAuto.csv";
	public final String VARIANCE_FILENAME = "resultsVariance.csv";
	private final String[] ROI_COLORS = {"#4dffff00","#3dff0000","#3d0000ff","#3d66ff66"};
	private final String[] DRAW_VELO_GROUP_CHOICES = {"Capillaries", "Arteries", "Veins", "Arteries/Veins", "All"};
	private final int[] DRAW_VELO_GROUP_CHOICES_VAL = {1, 2, 4, 6, 7};

	private String directory;
	private ImagePlus oriImp;
	private ImagePlus allStacksImp;
	private ImagePlus xtImp;
	private ImagePlus[] xtPreImp;
	private ImagePlus verifImp;
	private String[] xtPreImpMethod;
	private String oriROIsPath;
	private String oriName;
	private String xtPath;
	private String xtPathCrop;
	private String rtPath;
	private String rtManualPath;
	private String rtAutoPath;
	private String rtSlidingPath;
	private String rtSlidingAutoPath;
	private String varPath;
	private boolean showAll = false;
	
	private final int BORDER_SIZE = 2;
	private double[] varianceArray;
	private final int TABLE_PRECISION = 6;
	private ResultsTable rt;
	private ResultsTable var;
	private int[] rangeUID;
	private int idxUID;
	private int[] xyUID = new int[2]; // number of columns / lines
	private int[] whUID = new int[2]; // size of columns / lines
	
	private Object[] groupColorShortcutsCode = new Object[4];
	private String[] groupColorShortcutsKeys = {"1","2","3","4"};
	private String[] groupColorShortcutsName = {"RBCv SetGroup Capillaries", "RBCv SetGroup Arteries", "RBCv SetGroup Veins", "RBCv SetGroup Discard"};
	private Object[] setROIsShortcutsCode = new Object[4];
	private String[] setROIsShortcutsKeys = {"a","d","p","z"};
	private String[] setROIsShortcutsName = {"RBCv SetROIs Add Vertex", "RBCv SetROIs Delete Vertex", "RBCv SetROIs ToPolyline", "RBCv SetROIs Toggle Show All"};
	private Object[] verifShortcutsCode = new Object[4];
	private String[] verifShortcutsKeys = {"u","r","y","n"};
	private String[] verifShortcutsName = {"RBCv Verif Update Angle", "RBCv Verif Discard Angle", "RBCv Verif Validate ROI", "RBCv Verif Discard ROI"};
	private final String[] PRETREAT_METHODS = {"Demeaning","SobelX3","All"};
	//private final String[] pretreatMethods = {"Demeaning","SobelX2","SobelXY2","SobelX3","SobelXY3","SobelX5","SobelXY5","PrewittX3","PrewittXY3","All"};
	
	private final String FILE_COL_NAME = "File";
	private final String ROI_COL_NAME = "ROI";
	//private final String PREPROCESS_COL_NAME = "preprocess";
	private final String YPX_UM_COL_NAME = "yPx_um";
	private final String XPX_MS_COL_NAME = "xPx_ms";
	private final String LENGTH_COL_NAME = "length";
	private final String VELOCITY_COL_NAME = "velocity";
	private final String ANGLE_COL_NAME = "angle";
	private final String SEP_COL_NAME = "sep";
	private final String FWHM_COL_NAME = "FWHM";
	private final String CORREL_COL_NAME = "correl";
	private final String SD_EVOL_COL_NAME = "sd_evol";
	private final String GROUP_COL_NAME = "group";
	private final String VERIF_COL_NAME = "Verified";
	private final String VELOCITY_OK_COL_NAME = "velocity_ok";
	private final String ANGLE_OK_COL_NAME = "angle_ok";
	private final String RS_COEFF_COL_NAME = "rs_coeff";
	private final String VARIANCE_COL_NAME = "variance";
	private final String[] VERIF_METRICS = {SEP_COL_NAME,FWHM_COL_NAME,CORREL_COL_NAME,SD_EVOL_COL_NAME};
	private final boolean[] VERIF_METRICS_ABOVE = {true, false, true, true};
	
	protected RBCVProcess() {}
	
	public static RBCVProcess getInstance() {
		if (instance == null) {
			instance = new RBCVProcess();
		}
		return instance;
	}
	
	public void execSubFunction(String arg0) {
		if (arg0.matches("radon_transform")) {
			radonTransformDialog();
			return;
		}
		if (frame == null) {
			getFrame();
			return;
		}
		if (arg0.matches("setrois_add_vertex")) {
			addVertex();
		} else if (arg0.matches("setrois_delete_vertex")) {
			deleteVertex();
		} else if (arg0.matches("setrois_to_polyline")) {
			toPolyline();
		} else if (arg0.matches("setrois_toggle_show_all")) {
			toggleShowAll();
		} else if (arg0.matches("verif_update_angle")) {
			updateAngle();
		} else if (arg0.matches("verif_discard_angle")) {
			discardAngle();
		} else if (arg0.matches("verif_validate_roi")) {
			validateROI();
		} else if (arg0.matches("verif_unvalidate_roi")) {
			unvalidateROI();
		} else if (arg0.matches("setgroup_capillaries")) {
			setColorGroup(0);
		} else if (arg0.matches("setgroup_arteries")) {
			setColorGroup(1);
		} else if (arg0.matches("setgroup_veins")) {
			setColorGroup(2);
		} else if (arg0.matches("setgroup_discard")) {
			setColorGroup(3);
		} else {
			getFrame();
		}
	}
	
	public void getFrame() {
		if (frame == null) {
			if (depthValue.getItemCount() == 0) {
				for (int i = 0; i < 4; i++) {
					depthValue.add("Depth: "+i);
				}
			}
			if (!frameInitialized) {
				selectDirectoryBtn.addActionListener(this);
				depthValue.addItemListener(this);
				stackSuffixValue.setColumns(10);
				stackSuffixValue.addTextListener(this);
				stackSuffixValue.addActionListener(this);
				stackSuffixValue.addFocusListener(this);
				//shiftSuffixValue.setColumns(10);
				//shiftSuffixValue.addTextListener(this);
				//shiftSuffixValue.addFocusListener(this);
				frameInitialized = true;
				openStackBtn.addActionListener(this);
				saveROIsBtn.addActionListener(this);
				exportROIsBtn.addActionListener(this);
				xtGenBtn.addActionListener(this);
				measureVelocityBtn.addActionListener(this);
				measureVelocitySlidingBtn.addActionListener(this);
				manualVerifBtn.addActionListener(this);
				evaluateVerifBtn.addActionListener(this);
				automaticVerifBtn.addActionListener(this);
				rollingShutterBtn.addActionListener(this);
				velocityLimitBtn.addActionListener(this);
				drawVelocityBtn.addActionListener(this);
			}
			
			frame = new Frame("Red Blood Cell velocity v1.0.0");
			frame.setLayout(new GridBagLayout());
			
			Container stacksLayout = new Container();
			stacksLayout.setLayout(new GridBagLayout());
			addThingContainer(stacksLayout, Box.createVerticalStrut(6),		0, 0,	1, 101,	0, 0,	6, 0); // left margin
			addThingContainer(stacksLayout, Box.createVerticalStrut(6),		7, 0,	1, 101,	0, 0,	6, 0); // right margin
			addThingContainer(stacksLayout, Box.createHorizontalStrut(6),	1, 0,	2, 1,	0, 0,	0, 6); // top margin
			addThingContainer(stacksLayout, Box.createHorizontalStrut(3),	1, 100,	2, 1,	0, 0,	0, 3); // bottom margin
			addThingContainer(stacksLayout, selectDirectoryBtn,				1, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, directoryLabel,					1, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, directoryValue,					2, 2,	5, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, depthValue,						2, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, stacksListLabel,				1, 3,	1, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, stacksListValue,				2, 3,	5, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, stackSuffixLabel,				3, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(stacksLayout, stackSuffixValue,				4, 1,	1, 1,	1, 1,	0, 0);
			//addThingContainer(stacksLayout, shiftSuffixLabel,				5, 1,	1, 1,	1, 1,	0, 0);
			//addThingContainer(stacksLayout, shiftSuffixValue,				6, 1,	1, 1,	1, 1,	0, 0);
			
			Container roiLayout = new Container();
			roiLayout.setLayout(new GridBagLayout());
			addThingContainer(roiLayout, Box.createVerticalStrut(6),	0, 0,	1, 101,	0, 0,	6, 0); // left margin
			addThingContainer(roiLayout, Box.createVerticalStrut(6),	7, 0,	1, 101,	0, 0,	6, 0); // right margin
			addThingContainer(roiLayout, Box.createHorizontalStrut(3),	1, 100,	2, 1,	0, 0,	0, 3); // bottom margin
			addThingContainer(roiLayout, openStackBtn,					1, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, saveROIsBtn,					2, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, exportROIsBtn,					3, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, xtGenBtn,						1, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, measureVelocityBtn,			2, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, measureVelocitySlidingBtn,		3, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, manualVerifBtn,				1, 3,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, evaluateVerifBtn,				2, 3,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, automaticVerifBtn,				3, 3,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, velocityLimitBtn,				1, 4,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, rollingShutterBtn,				2, 4,	1, 1,	1, 1,	0, 0);
			addThingContainer(roiLayout, drawVelocityBtn,				3, 4,	1, 1,	1, 1,	0, 0);

			addThingContainer(frame, stacksLayout,			0, 0,	1, 1,	1, 1,	0, 0);
			addThingContainer(frame, roiLayout,				0, 1,	1, 1,	1, 1,	0, 0);
			
			//addThingContainer(frame, Box.createHorizontalStrut(6),	1, 85,	2, 1,	0, 0,	0, 6);
			
			//addThingContainer(frame, Box.createHorizontalStrut(6),	1, 101,	2, 1,	0, 0,	0, 6);
			
			frame.setLocation((int)framePosX, (int)framePosY+21);
	        frame.pack();
			//frame.setSize((int)frameSizeX, (int)frameSizeY);
	        frame.setVisible(true);
	        frame.addWindowListener(this);
	        
			Hashtable shortcuts = Menus.getShortcuts();
			for (int iKeys = 0; iKeys < setROIsShortcutsKeys.length; iKeys++) {
				int keyValue = Menus.convertShortcutToCode(setROIsShortcutsKeys[iKeys]);
				setROIsShortcutsCode[iKeys] = shortcuts.get(keyValue); // store shortcuts to put them back when verif is closed
				shortcuts.put(keyValue,setROIsShortcutsName[iKeys]);
			}
			for (int iKeys = 0; iKeys < groupColorShortcutsKeys.length; iKeys++) {
				int keyValue = Menus.convertShortcutToCode(groupColorShortcutsKeys[iKeys]);
				groupColorShortcutsCode[iKeys] = shortcuts.get(keyValue); // store shortcuts to put them back when verif is closed
				shortcuts.put(keyValue,groupColorShortcutsName[iKeys]);
			}
		} else {
			frame.toFront();
		}
		if (directory == null) {
			updateDirectory();
			if (directory == null) {
				frame.dispose();
				frame = null;
			}
		}
	}
	
	public void getVerif() {
		if (verif == null) {
			if (!verifInitialized) {
				verifInitialized = true;
				updateAngleBtn.addActionListener(this);
				discardAngleBtn.addActionListener(this);
				openStackVerifBtn.addActionListener(this);
				openXTRawVerifBtn.addActionListener(this);
				validateROIBtn.addActionListener(this);
				unvalidateROIBtn.addActionListener(this);
			}
			
			verif = new Frame("Verification");
			verif.setLayout(new GridBagLayout());
			
			addThingContainer(verif, Box.createVerticalStrut(6),		0, 0,	1, 101,	0, 0,	6, 0); // left margin
			addThingContainer(verif, Box.createVerticalStrut(6),		7, 0,	1, 101,	0, 0,	6, 0); // right margin
			addThingContainer(verif, Box.createHorizontalStrut(3),		1, 100,	2, 1,	0, 0,	0, 6); // bottom margin
			addThingContainer(verif, updateAngleBtn,					1, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, discardAngleBtn,					1, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, validateROIBtn,					2, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, unvalidateROIBtn,					2, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, openStackVerifBtn,					3, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, openXTRawVerifBtn,					3, 2,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, openVirtualCB,						4, 1,	1, 1,	1, 1,	0, 0);
			addThingContainer(verif, autoCloseCB,						4, 2,	1, 1,	1, 1,	0, 0);

			
			//addThingContainer(verif, Box.createHorizontalStrut(6),		1, 100,	2, 1,	0, 0,	0, 6);
			
			verif.setLocation((int)verifPosX, (int)verifPosY+21);
			verif.pack();
			//verif.setSize((int)frameSizeX, (int)frameSizeY);
			verif.setVisible(true);
			verif.addWindowListener(this);
			
			Hashtable shortcuts = Menus.getShortcuts();
			for (int iKeys = 0; iKeys < verifShortcutsKeys.length; iKeys++) {
				int keyValue = Menus.convertShortcutToCode(verifShortcutsKeys[iKeys]);
				verifShortcutsCode[iKeys] = shortcuts.get(keyValue); // store shortcuts to put them back when verif is closed
				shortcuts.put(keyValue,verifShortcutsName[iKeys]);
			}
		} else {
			verif.toFront();
		}
	}
	
	private void addThingContainer(
		Container f, Component b, 
		int gridx, int gridy, int gridwidth, int gridheight, 
		int weightx, int weighty, int ipadx,int ipady
	){
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.BOTH;
		c.gridx = gridx;
		c.gridy = gridy;
		c.gridwidth = gridwidth;
		c.gridheight = gridheight;
		c.weightx = weightx;
		c.weighty = weighty;
		c.ipadx = ipadx;
		c.ipady = ipady;
		f.add(b, c);
	}
	
	private void addVertex() {
		ImagePlus imp = IJ.getImage();
		Roi r = imp.getRoi();
		if (r == null) {
			new MessageDialog(new Frame(), "Incorrect selection", "No selection to add vertex to.");
			return;
		}
		int selType = r.getType();
		if (selType == Roi.LINE || selType == Roi.POLYLINE || selType == Roi.FREELINE) {
			Point cursorPos = imp.getCanvas().getCursorLoc();
			float[] xPoints = r.getFloatPolygon().xpoints;
			float[] yPoints = r.getFloatPolygon().ypoints;
			if (selType == Roi.LINE && xPoints.length > 2) { // if selType == 5, xPoints has length 10 instead of 2, paded with 0s...
				float[] tmp = {xPoints[0], xPoints[1], yPoints[0], yPoints[1]};
				xPoints = new float[2];
				yPoints = new float[2];
				xPoints[0] = tmp[0];
				xPoints[1] = tmp[1];
				yPoints[0] = tmp[2];
				yPoints[1] = tmp[3];
			}
			int roiID = getRM().getSelectedIndex();
			Color color = r.getStrokeColor();
			float width = r.getStrokeWidth();
			double mindist = Double.POSITIVE_INFINITY;
			int index = -1;
			int last = xPoints.length-1;
			for (int i = 0; i < xPoints.length; i++) {
				double dist = Math.sqrt(Math.pow(xPoints[i]-cursorPos.x,2)+Math.pow(yPoints[i]-cursorPos.y,2));
				if (mindist > dist) {
					mindist = dist;
					index = i;
				}
			}
			float V0x = cursorPos.x-xPoints[index];
			float V0y = cursorPos.y-yPoints[index];
			if (index == 0) {
				if (Math.sqrt(Math.pow(xPoints[0]-xPoints[1],2)+Math.pow(yPoints[0]-yPoints[1],2)) > 
						Math.sqrt(Math.pow(xPoints[1]-cursorPos.x,2)+Math.pow(yPoints[1]-cursorPos.y,2))) {
					index++;
				}
			} else if (index == last) {
				if (Math.sqrt(Math.pow(xPoints[last]-xPoints[last-1],2)+Math.pow(yPoints[last]-yPoints[last-1],2)) < 
						Math.sqrt(Math.pow(xPoints[last-1]-cursorPos.x,2)+Math.pow(yPoints[last-1]-cursorPos.y,2))) {
					index++;
				}
			} else {
				float V1x = xPoints[index-1] - xPoints[index];
				float V1y = yPoints[index-1] - yPoints[index];
				float V2x = xPoints[index+1] - xPoints[index];
				float V2y = yPoints[index+1] - yPoints[index];
				// theta = arccos ( V1*V2 / ( ||V1||*||V2|| ) ) where V1 and V2 are vectors
				double theta1 = Math.acos((V1x*V0x+V1y*V0y) / (Math.sqrt(V1x*V1x+V1y*V1y)*Math.sqrt(V0x*V0x+V0y*V0y)));
				double theta2 = Math.acos((V2x*V0x+V2y*V0y) / (Math.sqrt(V2x*V2x+V2y*V2y)*Math.sqrt(V0x*V0x+V0y*V0y)));
				if (Math.abs(theta1) > Math.abs(theta2)) {
					index++;
				}
			}
			float[] xCoord = new float[xPoints.length+1];
			float[] yCoord = new float[yPoints.length+1];
			System.arraycopy(xPoints, 0, xCoord, 0, index);
			System.arraycopy(yPoints, 0, yCoord, 0, index);
			System.arraycopy(xPoints, index, xCoord, index+1, xPoints.length-index);
			System.arraycopy(yPoints, index, yCoord, index+1, yPoints.length-index);
			xCoord[index] = cursorPos.x;
			yCoord[index] = cursorPos.y;
			r = new PolygonRoi(xCoord, yCoord, Roi.POLYLINE);
			r.setStrokeColor(color);
			r.setStrokeWidth(width);
			imp.setRoi(r);
			if (roiID >= 0) {
				getRM().runCommand(imp,"Update");
				getRM().select(roiID);
			}
		} else {
			new MessageDialog(new Frame(), "Incorrect selection", "Could not add vertex, selection is not Straight, Segmented or Freehand Line.");
		}
	}
	
	private void deleteVertex() {
		ImagePlus imp = IJ.getImage();
		Roi r = imp.getRoi();
		if (r == null) {
			new MessageDialog(new Frame(), "Incorrect selection", "No selection to remove vertex from.");
			return;
		}
		int selType = r.getType();
		if (selType == Roi.POLYLINE || selType == Roi.FREELINE) {
			Point cursorPos = imp.getCanvas().getCursorLoc();
			float[] xPoints = r.getFloatPolygon().xpoints;
			float[] yPoints = r.getFloatPolygon().ypoints;
			Color color = r.getStrokeColor();
			float width = r.getStrokeWidth();
			if (xPoints.length > 2) {
				int roiID = getRM().getSelectedIndex();
				double mindist = Double.POSITIVE_INFINITY;
				int index = -1;
				for (int i = 0; i < xPoints.length; i++) {
					double dist = Math.sqrt(Math.pow(xPoints[i]-cursorPos.x,2)+Math.pow(yPoints[i]-cursorPos.y,2));
					if (mindist > dist) {
						mindist = dist;
						index = i;
					}
				}
				float[] xCoord = new float[xPoints.length-1];
				float[] yCoord = new float[yPoints.length-1];
				System.arraycopy(xPoints, 0, xCoord, 0, index);
				System.arraycopy(yPoints, 0, yCoord, 0, index);
				System.arraycopy(xPoints, index+1, xCoord, index, xPoints.length-index-1);
				System.arraycopy(yPoints, index+1, yCoord, index, yPoints.length-index-1);
				r = new PolygonRoi(xCoord, yCoord, Roi.POLYLINE);
				r.setStrokeColor(color);
				r.setStrokeWidth(width);
				imp.setRoi(r);
				if (roiID >= 0) {
					getRM().runCommand(imp,"Update");
					getRM().select(roiID);
				}
			} else {
				new MessageDialog(new Frame(), "Incorrect selection", "Could not remove vertex, selection only has two of them.");
			}
		} else {
			new MessageDialog(new Frame(), "Incorrect selection", "Could not remove vertex, selection is not Segmented or Freehand Line.");
		}
		
	}
	
	private void toPolyline() {
		toPolyline(20);
	}
	
	private void toPolyline(double spacing) {
		toPolyline(spacing, true, true);
	}
	
	private void toPolyline(double spacing, boolean verbose, boolean updateRM) {
		toPolyline(spacing, verbose, updateRM, IJ.getImage());
	}
	
	private void toPolyline(double spacing, boolean verbose, boolean updateRM, ImagePlus imp) {
		// change current ROI to a polyline with equal segments of size approximately "spacing"
		Roi r = imp.getRoi();
		if (r == null) {
			if (verbose) new MessageDialog(new Frame(), "Selection needed", "No selection to convert to polyline.");
			return;
		}
		int selType = r.getType();
		if (selType == Roi.LINE || selType == Roi.POLYLINE || selType == Roi.FREELINE) {
			float[] xIn = r.getFloatPolygon().xpoints;
			float[] yIn = r.getFloatPolygon().ypoints;
			if (selType == Roi.LINE && xIn.length > 2) { // if selType == 5, xPoints has length 10 instead of 2, paded with 0s...
				float[] tmp = {xIn[0], xIn[1], yIn[0], yIn[1]};
				xIn = new float[2];
				xIn = new float[2];
				xIn[0] = tmp[0];
				xIn[1] = tmp[1];
				yIn[0] = tmp[2];
				yIn[1] = tmp[3];
			}
			Color color = r.getStrokeColor();
			float width = r.getStrokeWidth();
			int roiID = getRM().getSelectedIndex();
			int nIn = xIn.length;
			double[] xyIn = new double[nIn];
			for (int iIn = 1; iIn < nIn; iIn++) {
				xyIn[iIn] = xyIn[iIn-1] + Math.sqrt(Math.pow(xIn[iIn]-xIn[iIn-1],2)+Math.pow(yIn[iIn]-yIn[iIn-1],2));
			}
			
			int nOut = (int)Math.round(xyIn[nIn-1]/spacing);
			spacing = xyIn[nIn-1]/nOut;
			float[] xOut = new float[nOut+1];
			float[] yOut = new float[nOut+1];
			xOut[0] = xIn[0];
			yOut[0] = yIn[0];
			
			int idx = 1;
			for (int i = 0; i < nOut; i++) {
				double pos = i*spacing;
				while (xyIn[idx] < pos) {
					idx++;
				}
				double coeff1 = (xyIn[idx]-pos) / (xyIn[idx]-xyIn[idx-1]);
				xOut[i] = (float)(coeff1 * xIn[idx-1] + (1-coeff1) * xIn[idx]);
				yOut[i] = (float)(coeff1 * yIn[idx-1] + (1-coeff1) * yIn[idx]);
			}
			xOut[nOut] = xIn[nIn-1];
			yOut[nOut] = yIn[nIn-1];
			r = new PolygonRoi(xOut, yOut, Roi.POLYLINE);
			r.setStrokeColor(color);
			r.setStrokeWidth(width);
			imp.setRoi(r);
			if (roiID >= 0 && updateRM) {
				getRM().runCommand(imp,"Update");
				getRM().select(roiID);
			}
		}
	}
		
	private void toggleShowAll() {
		showAll = !showAll;
		if (showAll) {
			getRM().runCommand("Show All with labels");
		} else {
			getRM().runCommand("Show None");
		}
	}
	
	private void updateDirectory() {
		DirectoryChooser dc = new DirectoryChooser("Choose directory containing image stacks");
		if (dc.getDirectory() != null) {
			if (directory != null && dc.getDirectory().matches(directory)) return;
			directory = dc.getDirectory();
			rtPath = directory+VELOCITY_FILENAME;
			rtManualPath = directory+VELOCITY_MANUAL_FILENAME;
			rtAutoPath = directory+VELOCITY_AUTO_FILENAME;
			rtSlidingPath = directory+VELOCITY_SLIDING_FILENAME;
			rtSlidingAutoPath = directory+VELOCITY_SLIDING_AUTO_FILENAME;
			varPath = directory+VARIANCE_FILENAME;
			directoryValue.setText(directory);
			updateStacksList();
			if (oriImp != null) {
				oriImp.close();
				oriImp = null;
			}
			if (allStacksImp != null) {
				allStacksImp.close();
				allStacksImp = null;
			}
		}
	}
	
	private void updateStacksList() {
		if (directory != null) {
			final File directoryFile = new File(directory);
			final int depth = depthValue.getSelectedIndex();
			ArrayList<String> fileList = getFilesDepthEnding(directoryFile, depth, stackSuffixValue.getText());
			stacksListValue.removeAll();
			for (String file : fileList) {
				if (file.endsWith(stackSuffixValue.getText())) {
					stacksListValue.add(file);
				}
			}
		}
	}
	
	public static ArrayList<String> getFilesDepthEnding(final File directory, final int depth, final String ending) {
		ArrayList<String> res = new ArrayList<String>();
	    for (final File fileEntry : directory.listFiles()) {
	    	boolean isDir = fileEntry.isDirectory();
	        if (isDir && depth > 0) {
	            ArrayList<String> tempRes = getFilesDepthEnding(fileEntry, depth-1, ending);
	            for (String temp : tempRes) {
	            	res.add(fileEntry.getName()+File.separator+temp);
	            }
	        } else if (!isDir && depth == 0) {
        		res.add(fileEntry.getName());
	        }
	    }
		return res;
	}
	
	private void openStack(boolean openVirtual, boolean openOneSlicePerStack) {
		if (!new File(directory+stacksListValue.getSelectedItem()).exists()) {
			new MessageDialog(new Frame(), "Missing stack", "No stack to open.");
			return;
		}
		String filePath = directory+stacksListValue.getSelectedItem();
		String fileFormat = Opener.getFileFormat(filePath);
		boolean openClassic = !fileFormat.matches("unknown") && !fileFormat.matches("txt");
		if (openVirtual) {
			if (openClassic) {
				oriImp = IJ.openVirtual(filePath);
			} else {
				IJ.run("Bio-Formats Importer", "open=["+filePath+
						"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
				oriImp = IJ.getImage();
			}
		} else {
			if (openClassic) {
				oriImp = IJ.openImage(filePath);
			} else {
				IJ.run("Bio-Formats Importer", "open=["+filePath+
						"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
				oriImp = IJ.getImage();
			}
		}
		oriImp.show();
		if (openOneSlicePerStack && (allStacksImp == null || !allStacksImp.isVisible())) {
			allStacksImp = new ImagePlus();
			for (int i = 0; i < stacksListValue.getItemCount(); i++) {
				ImagePlus tmpImp;
				if (openClassic) {
					tmpImp = IJ.openVirtual(directory+stacksListValue.getItem(i));
				} else {
					IJ.run("Bio-Formats Importer", "open=["+directory+stacksListValue.getItem(i)+
							"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
					tmpImp = IJ.getImage();
				}
				int middleSlice = (int)Math.round(tmpImp.getNSlices()/2.0);
				tmpImp.setSlice(middleSlice);
				if (i == 0) {
					allStacksImp = tmpImp.crop();
				} else {
					ImagePlus tmpImp2 = tmpImp.crop();
					allStacksImp = Concatenator.run(allStacksImp, tmpImp2);
				}
				tmpImp.close();
			}
			allStacksImp.setTitle("AllStacks");
			allStacksImp.show();
		}
		setPaths(directory+stacksListValue.getSelectedItem());
		openImageROIs(true, true);
	}
	
	private void saveROIs() {
		RoiManager rm = getRM();
		if (rm.getCount() > 0) {
			rm.deselect();
			File f = new File(oriROIsPath);
			f = f.getParentFile();
			if (!f.exists()) {
				f.mkdir();
			}
			rm.runCommand("Save", oriROIsPath);
		}
	}
	
	private void exportROIsToOtherStacks() {
		if (emptyListAbort()) return;
		if (getRM().getCount() == 0) {
			new MessageDialog(new Frame(), "Missing ROI", "No ROI to export. Aborting.");
			return;
		}
		YesNoCancelDialog d = new YesNoCancelDialog(new Frame(), "", "Are you sure you want to export ROIs to all other\n"+
				"stacks, possibly overwritting exising ones?");
		if (d.yesPressed()) {
			RoiManager rm = getRM();
			if (rm.getCount() > 0) {
				rm.deselect();
				for (int i = 0; i < stacksListValue.getItemCount(); i++) {
					setPaths(directory+stacksListValue.getItem(i));
					rm.runCommand("Save", oriROIsPath);
				}
			}
		}
	}
	
	private void setPaths(String imagePath) {
		File img = new File(imagePath);
		String parentPath = img.getParent();
		String imgName = img.getName();
		int lastIdx = imgName.lastIndexOf(".");
		oriName = imgName.substring(0, lastIdx);
		oriROIsPath = parentPath+File.separator+"ROISets"+File.separator+oriName+"_ROIs.zip";
		File f = new File(oriROIsPath);
		f = f.getParentFile();
		if (!f.exists()) {
			f.mkdir();
		}
		xtPath = parentPath+File.separator+"XT"+File.separator;
		File xtFile = new File(xtPath);
		if (!xtFile.exists()) {
			xtFile.mkdir();
		}
		xtPathCrop = parentPath+File.separator+"XTcropParms"+File.separator;
		File xtFileCrop = new File(xtPathCrop);
		if (!xtFileCrop.exists()) {
			xtFileCrop.mkdir();
		}
	}
	
	private boolean emptyListAbort() {
		if (stacksListValue.getItemCount() < 1) {
			new MessageDialog(new Frame(), "Empty list", "Stack list is empty. Aborting.");
			return true;
		} else {
			return false;
		}
	}
	
	private void xtGen() {
		if (emptyListAbort()) return;
		GenericDialog gd = new GenericDialog("XT Gen");
		gd.addNumericField("Gaussian sigma", Prefs.get("RBCv.gaussSigma", 0), 1, 8, "px");
		gd.addNumericField("Background width factor", Prefs.get("RBCv.backgroundWidthFactor", 2));
		gd.showDialog();
		if (gd.wasCanceled()) return;
		double gaussSigma = gd.getNextNumber();
		Prefs.set("RBCv.gaussSigma", gaussSigma);
		double backgroundWidthFactor = gd.getNextNumber();
		Prefs.set("RBCv.backgroundWidthFactor", backgroundWidthFactor);
		generateXTImages(gaussSigma, backgroundWidthFactor);
	}
	
	private void generateXTImages(double gaussSigma, double backgroundWidthFactor) {
		if (emptyListAbort()) return;
		boolean removeBackground = backgroundWidthFactor > 1;
		long startTime = System.currentTimeMillis();
		if (oriImp != null) oriImp.close();
		RoiManager rm = getRM();
		rm.reset();
		int nImages = stacksListValue.getItemCount();
		for (int iImage = 0; iImage < nImages; iImage++) {
			setPaths(directory+stacksListValue.getItem(iImage));
			openImageROIs(true, true);
			int nROI = rm.getCount();
			IJ.log("Image "+(iImage+1)+"/"+nImages+", "+nROI+" ROIs. "+((System.currentTimeMillis()-startTime)/1000.0)+"s elapsed.");
			
			oriImp = IJ.openImage(directory+stacksListValue.getItem(iImage));
			if (gaussSigma > 0) {
				IJ.run(oriImp, "Gaussian Blur...", "sigma="+gaussSigma+" stack");
			}
			for (int iROI = 0; iROI < nROI; iROI++) {
				rm.select(oriImp, iROI);
				//if (oriImp.getRoi().getType() == Roi.FREELINE) {
				//	toPolyline(5, false, false, oriImp);
				//}
				ProfilePlot profile = new ProfilePlot(oriImp);
				double[] values = profile.getProfile();
				int size = values.length;
				int slices = oriImp.getNSlices() * oriImp.getNFrames();
				float[] data = new float[size*slices];
				for (int iSlice = 0; iSlice < slices; iSlice++) {
					oriImp.setSliceWithoutUpdate(iSlice+1);
					profile = new ProfilePlot(oriImp);
					values = profile.getProfile();
					for (int iSize = 0; iSize < size; iSize++) {
						data[iSize*slices+iSlice] = (float)values[iSize];
					}
				}
				if (removeBackground) {
					Roi r = oriImp.getRoi();
					float widthOri = r.getStrokeWidth();
					double widthNew = Math.max(widthOri+1, Math.round(backgroundWidthFactor * widthOri));
					float widthRatio = (float)widthNew / widthOri;
					r.setStrokeWidth(widthNew);
					oriImp.setRoi(r);
					for (int iSlice = 0; iSlice < slices; iSlice++) {
						oriImp.setSliceWithoutUpdate(iSlice+1);
						profile = new ProfilePlot(oriImp);
						values = profile.getProfile();
						for (int iSize = 0; iSize < size; iSize++) {
							data[iSize*slices+iSlice] = data[iSize*slices+iSlice] - (data[iSize*slices+iSlice]-(float)values[iSize]*widthRatio)/(widthRatio-1);
						}
					}
				}
				
				ImagePlus xtImp = IJ.createImage(ROI_COL_NAME+(iROI+1)+"init", "32-bit", slices, size, 1);
				xtImp.getProcessor().setPixels(data);
				Rectangle r = getNotNullBoundingBox(xtImp);
				if (r.y != 0 || r.height != xtImp.getHeight()) {
					if (r.height == 0) r.height = 1;
					xtImp.setRoi(0, r.y, xtImp.getWidth(), r.height);
					xtImp = xtImp.crop();
					String pathXtCropFile = xtPathCrop+"XTcrop"+(iROI+1)+"_"+oriName+".txt";
					String textOut = ""+r.y+"-"+(r.y+r.height-1);
					IJ.saveString(textOut, pathXtCropFile);
				}
				xtImp.getProcessor().resetMinAndMax();
				//xtImp.setProcessor(xtImp.getProcessor().convertToShortProcessor());
				IJ.saveAs(xtImp, "tif", xtPath+"XT"+(iROI+1)+"_"+oriName);
				xtImp.close();
			}
			oriImp.close();
		}
		IJ.log("XT images generation done. "+((System.currentTimeMillis()-startTime)/1000.0)+"s elapsed.");
	}
	
	private void measureVelocity(boolean slidingWindow) {
		if (emptyListAbort()) return;
		rt = new ResultsTable();
		rt.setPrecision(TABLE_PRECISION);
		rt.showRowNumbers(false);
		var = new ResultsTable();
		var.setPrecision(TABLE_PRECISION);
		var.showRowNumbers(false);
		int rtIdx = 0;
		double slidingWindowTimeDuration = Prefs.get("RBCv.slidingWindowTimeDuration", 0.50);
		double slidingWindowTimeStep = Prefs.get("RBCv.slidingWindowTimeStep", 0.10);
		int slidingWindowPxDuration = 0;
		int slidingWindowPxStep = 0;
		
		GenericDialog gd = new GenericDialog("Measure velocity");
		gd.addMessage("Global parameters");
		gd.addNumericField("Pixel size", Prefs.get("RBCv.pixel_size_um", 1.3), 3, 8, "um");
		gd.addNumericField("Time between images", Prefs.get("RBCv.time_between_images_ms", 10), 3, 8, "ms");
		gd.addMessage("Pretreat parameters:");
		gd.addNumericField("Gaussian sigma x", Prefs.get("RBCv.gaussSigmaX", 0), 1, 8, "px");
		//gd.addChoice("Pretreat method: ", pretreatMethods, Prefs.get("RBCv.pretreatMethod", "Demeaning")); // SobelX3 only to simplify
		if (slidingWindow) {
			gd.addMessage("Sliding window parameters:");
			gd.addNumericField("Window size", slidingWindowTimeDuration, 2, 8, "s");
			gd.addNumericField("Window step", slidingWindowTimeStep, 2, 8, "s");
		} else {
			gd.addCheckbox("Preview pretreat only", false);
		}
		gd.showDialog();
		if (gd.wasCanceled()) return;
		boolean previewOnly = false;
		double pixel_size_um = gd.getNextNumber();
		Prefs.set("RBCv.pixel_size_um", pixel_size_um);
		double time_between_images_ms = gd.getNextNumber();
		Prefs.set("RBCv.time_between_images_ms", time_between_images_ms);
		double gaussSigmaX = gd.getNextNumber();
		Prefs.set("RBCv.gaussSigmaX", gaussSigmaX);
		String pretreatMethod = "SobelX3"; // SobelX3 only to simplify
		//String pretreatMethod = gd.getNextChoice();
		//Prefs.set("RBCv.pretreatMethod", pretreatMethod);
		if (slidingWindow) {
			slidingWindowTimeDuration = gd.getNextNumber();
			Prefs.set("RBCv.slidingWindowTimeDuration", slidingWindowTimeDuration);
			slidingWindowTimeStep = gd.getNextNumber();
			Prefs.set("RBCv.slidingWindowTimeStep", slidingWindowTimeStep);
			slidingWindowPxDuration = (int)Math.ceil(slidingWindowTimeDuration / (time_between_images_ms / 1000));
			slidingWindowPxStep = (int)Math.ceil(slidingWindowTimeStep / (time_between_images_ms / 1000));
		} else {
			previewOnly = gd.getNextBoolean();
		}
		
		if (pretreatMethod == "All") {
			xtPreImp = new ImagePlus[PRETREAT_METHODS.length-1];
			xtPreImpMethod = new String[PRETREAT_METHODS.length-1];
		} else {
			xtPreImp = new ImagePlus[1];
			xtPreImpMethod = new String[1];
		}
		
		final double[] radonAngles = new double[180];
		double[] radonAngles2 = new double[201];
		for (int iStep = 0; iStep < 180; iStep++) {
			radonAngles[iStep] = iStep;
		}
		long startTime = System.currentTimeMillis();
		int nImages = stacksListValue.getItemCount();
		ResultsTable varTmp = new ResultsTable();
		for (int iImage = 0; iImage < nImages; iImage++) {
			setPaths(directory+stacksListValue.getItem(iImage));
			openImageROIs(true, true);
			int nROI = getRM().getCount();
			IJ.log("Image "+(iImage+1)+"/"+nImages+", "+nROI+" ROIs. "+((System.currentTimeMillis()-startTime)/1000.0)+"s elapsed.");
			for (int iROI = 0; iROI < nROI; iROI++) {
				String xtFullPath = xtPath+"XT"+(iROI+1)+"_"+oriName+".tif";
				if (new File(xtFullPath).exists()) {
					int roiColorGroup = getColorGroup(getRM().getRoi(iROI));
					xtImp = IJ.openImage(xtFullPath);
					if (xtImp.getType() != ImagePlus.GRAY32) {
						xtImp.setProcessor(xtImp.getProcessor().convertToFloat());
					}
					if (gaussSigmaX > 0) IJ.run(xtImp, "Gaussian Blur 3D...", "x="+gaussSigmaX+" y=0 z=0");
					int idx = 0;
					if (pretreatMethod == "All") {
						for (int iMethod = 0; iMethod < PRETREAT_METHODS.length-1; iMethod++) {
							xtPreImpMethod[idx] = PRETREAT_METHODS[iMethod];
							xtPreImp[idx++] = pretreatXTImage(xtImp, PRETREAT_METHODS[iMethod], false);
						}
					} else {
						xtPreImpMethod[idx] = pretreatMethod;
						xtPreImp[idx++] = pretreatXTImage(xtImp, pretreatMethod, false);
					}
					
					if (previewOnly) {
						xtImp.show();
						for (int iMethod = 0; iMethod < idx; iMethod++) {
							xtPreImp[iMethod].show();
						}
						return;
					}
					
					for (int iMethod = 0; iMethod < idx; iMethod++) {
						if (slidingWindow) {
							int nWindow = (xtImp.getWidth() - slidingWindowPxDuration) / slidingWindowPxStep + 1;
							for (int iWindow = 0; iWindow < nWindow; iWindow++) {
								xtPreImp[iMethod].setRoi(new Rectangle(iWindow*slidingWindowPxStep,0,slidingWindowPxDuration,xtImp.getHeight()));
								ImagePlus cropImp = xtPreImp[iMethod].crop();
								
								double velocity = Double.NaN;
								double centerVariance = Double.NaN;
								double separability = Double.NaN;
								double fwhm = Double.NaN;
								double correlation = Double.NaN;
								double[] sdMeasures = {Double.NaN, Double.NaN};
								if (roiColorGroup != 3) {
									ImageStatistics is;
									is = cropImp.getStatistics(ImageStatistics.STD_DEV);
									sdMeasures[0] = is.stdDev;
									if (xtPreImp[iMethod].getHeight() > 1) {
										ImagePlus tmpMeasure = cropImp.resize(cropImp.getWidth()/2, cropImp.getHeight()/2, "average");
										is = tmpMeasure.getStatistics(ImageStatistics.STD_DEV);
										sdMeasures[1] = is.stdDev;
										tmpMeasure.close();
									} else {
										sdMeasures[1] = sdMeasures[0];
									}
									
									ImagePlus radonImp = radonTransform(cropImp, radonAngles);
									double maxVar = computeNormalizedVerticalVariance(radonImp, -1);
									centerVariance = secondOrderPolynomialMax(radonAngles,varianceArray);
									separability = getSeparability(varianceArray);
									fwhm = getFullWidthHalfMax(varianceArray);
	
									if (!Double.isNaN(centerVariance)) {
										for (int iAngle = 0; iAngle < radonAngles2.length; iAngle++) {
											radonAngles2[iAngle] = centerVariance-1+iAngle/100.0;
										}
										radonImp = radonTransform(cropImp, radonAngles2);
										computeNormalizedVerticalVariance(radonImp, maxVar);
										centerVariance = secondOrderPolynomialMax(radonAngles2,varianceArray);
									}
									velocity = pixel_size_um*1000/(time_between_images_ms*Math.tan(centerVariance*Math.PI/180));
									correlation = getCorrelation(cropImp, centerVariance);
									radonImp.close();
								}
								rt.setValue(FILE_COL_NAME, rtIdx, stacksListValue.getItem(iImage));
								rt.setValue(ROI_COL_NAME, rtIdx, iROI+1);
								//rt.setValue(PREPROCESS_COL_NAME, rtIdx, xtPreImpMethod[iMethod]);
								rt.setValue("t", rtIdx, (slidingWindowPxDuration/2.0+iWindow*slidingWindowPxStep)*time_between_images_ms/1000.0);
								rt.setValue(YPX_UM_COL_NAME, rtIdx, pixel_size_um);
								rt.setValue(XPX_MS_COL_NAME, rtIdx, time_between_images_ms);
								rt.setValue(LENGTH_COL_NAME, rtIdx, cropImp.getHeight());
								rt.setValue(VELOCITY_COL_NAME, rtIdx, velocity);
								rt.setValue(ANGLE_COL_NAME, rtIdx, centerVariance);
								rt.setValue(SEP_COL_NAME, rtIdx, separability);
								rt.setValue(FWHM_COL_NAME, rtIdx, fwhm);
								rt.setValue(CORREL_COL_NAME, rtIdx, correlation);
								rt.setValue(SD_EVOL_COL_NAME, rtIdx, sdMeasures[1]/sdMeasures[0]);
								rt.setValue(GROUP_COL_NAME, rtIdx, roiColorGroup);
								rtIdx++;
								cropImp.close();
							}
						} else {
							double velocity = Double.NaN;
							double centerVariance = Double.NaN;
							double separability = Double.NaN;
							double fwhm = Double.NaN;
							double correlation = Double.NaN;
							double[] sdMeasures = {Double.NaN, Double.NaN};
							if (roiColorGroup != 3) {
								ImageStatistics is = xtPreImp[iMethod].getStatistics(ImageStatistics.STD_DEV);
								sdMeasures[0] = is.stdDev;
								if (xtPreImp[iMethod].getHeight() > 1) {
									ImagePlus tmpMeasure = xtPreImp[iMethod].resize(xtPreImp[iMethod].getWidth()/2, xtPreImp[iMethod].getHeight()/2, "average");
									is = tmpMeasure.getStatistics(ImageStatistics.STD_DEV);
									sdMeasures[1] = is.stdDev;
									tmpMeasure.close();
								} else {
									sdMeasures[1] = sdMeasures[0];
								}
								ImagePlus radonImp = radonTransform(xtPreImp[iMethod], radonAngles);
								double maxVar = computeNormalizedVerticalVariance(radonImp, -1);
								centerVariance = secondOrderPolynomialMax(radonAngles,varianceArray);
								varTmp.reset();
								for (int iAngle = 0; iAngle < radonAngles.length; iAngle++) {
									int newLine = varTmp.size();
									varTmp.setValue(ANGLE_COL_NAME, newLine, radonAngles[iAngle]);
									varTmp.setValue(VARIANCE_COL_NAME, newLine, varianceArray[iAngle]);
								}
								separability = getSeparability(varianceArray);
								fwhm = getFullWidthHalfMax(varianceArray);
								
								if (!Double.isNaN(centerVariance)) {
									for (int iAngle = 0; iAngle < radonAngles2.length; iAngle++) {
										radonAngles2[iAngle] = centerVariance-1+iAngle/100.0;
									}
									radonImp = radonTransform(xtPreImp[iMethod], radonAngles2);
									computeNormalizedVerticalVariance(radonImp, maxVar);
									centerVariance = secondOrderPolynomialMax(radonAngles2,varianceArray);
									for (int iAngle = 0; iAngle < radonAngles2.length; iAngle++) {
										int newLine = varTmp.size();
										varTmp.setValue(ANGLE_COL_NAME, newLine, radonAngles2[iAngle]);
										varTmp.setValue(VARIANCE_COL_NAME, newLine, varianceArray[iAngle]);
									}
									varTmp.sort(ANGLE_COL_NAME);
								}
								velocity = pixel_size_um*1000/(time_between_images_ms*Math.tan(centerVariance*Math.PI/180));
								correlation = getCorrelation(xtPreImp[iMethod], centerVariance);
								radonImp.close();
							}
							String angleHeader = "angle_"+oriName+"_ROI"+(iROI+1)+"_"+xtPreImpMethod[iMethod];
							String varHeader = "variance_"+oriName+"_ROI"+(iROI+1)+"_"+xtPreImpMethod[iMethod];
							for (int iAngle = 0; iAngle < varTmp.size(); iAngle++) {
								var.setValue(angleHeader, iAngle, varTmp.getValue(ANGLE_COL_NAME,iAngle));
								var.setValue(varHeader, iAngle, varTmp.getValue(VARIANCE_COL_NAME,iAngle));
							}
							
							rt.setValue(FILE_COL_NAME, rtIdx, stacksListValue.getItem(iImage));
							rt.setValue(ROI_COL_NAME, rtIdx, iROI+1);
							//rt.setValue(PREPROCESS_COL_NAME, rtIdx, xtPreImpMethod[iMethod]);
							rt.setValue(YPX_UM_COL_NAME, rtIdx, pixel_size_um);
							rt.setValue(XPX_MS_COL_NAME, rtIdx, time_between_images_ms);
							rt.setValue(LENGTH_COL_NAME, rtIdx, xtPreImp[iMethod].getHeight());
							rt.setValue(VELOCITY_COL_NAME, rtIdx, velocity);
							rt.setValue(ANGLE_COL_NAME, rtIdx, centerVariance);
							rt.setValue(SEP_COL_NAME, rtIdx, separability);
							rt.setValue(FWHM_COL_NAME, rtIdx, fwhm);
							rt.setValue(CORREL_COL_NAME, rtIdx, correlation);
							rt.setValue(SD_EVOL_COL_NAME, rtIdx, sdMeasures[1]/sdMeasures[0]);
							rt.setValue(GROUP_COL_NAME, rtIdx, roiColorGroup);
							//rt.setValue(VERIF_COL_NAME, rtIdx, "0");
							rtIdx++;
						}
					}
				}
			}
		}
		if (rt.getHeadings().length < 1) {
			new MessageDialog(new Frame(), "Empty table", "No column in tables. Execute XT Gen first.");
			return;
		}
		var.save(varPath);
		//rt.sort(PREPROCESS_COL_NAME);
		rt.sort(ROI_COL_NAME);
		if (slidingWindow) {
			rt.save(rtSlidingPath);
			rt.show(VELOCITY_SLIDING_FILENAME);
		} else {
			rt.save(rtPath);
			rt.show(VELOCITY_FILENAME);
		}
		IJ.log("Velocity measures done. "+((System.currentTimeMillis()-startTime)/1000.0)+"s elapsed.");
	}
	
	public static double getSeparability(double[] varianceArray) {
		// separability = variance_max / mean(variance)
		double sumVariance = 0;
		double maxValue = Double.NEGATIVE_INFINITY;
		for (int iVariance = 0; iVariance < varianceArray.length; iVariance++) {
			sumVariance += varianceArray[iVariance];
			if (varianceArray[iVariance] > maxValue) {
				maxValue = varianceArray[iVariance];
			}
		}
		return maxValue*varianceArray.length/sumVariance;
	}
	
	public static double getFullWidthHalfMax(double[] array) {
		// FWHM in terms in number of array elements
		double maxValue = Double.NEGATIVE_INFINITY;
		int maxIndex = 0;
		for (int i = 0; i < array.length; i++) {
			if (array[i] > maxValue) {
				maxValue = array[i];
				maxIndex = i;
			}
		}
		double halfMaxValue = maxValue / 2.0;
		double[] permut = new double[array.length+1];
		System.arraycopy(array, 0, permut, array.length-maxIndex, maxIndex+1);
		System.arraycopy(array, maxIndex, permut, 0, array.length-maxIndex);
		double fwhm = 0;
		for (int i = 1; i < array.length; i++) {
			if (permut[i] < halfMaxValue) {
				fwhm += i-1+(permut[i-1]-halfMaxValue)/(permut[i-1]-permut[i]);
				break;
			}
		}
		for (int i = array.length-1; i > 0; i--) {
			if (permut[i] < halfMaxValue) {
				fwhm += array.length-1-i+(permut[i+1]-halfMaxValue)/(permut[i+1]-permut[i]);
				break;
			}
		}
		if (fwhm == 0) {
			fwhm = array.length;
		}
		return fwhm;
	}
	
	private double getCorrelation(ImagePlus cropImp, double centerVariance) {
		ImagePlus rotImp = cropImp.duplicate();
		IJ.run(rotImp, "Rotate... ", "angle="+centerVariance+" grid=1 interpolation=Bicubic enlarge");
		int widthOut = rotImp.getWidth();
		int heightOut = rotImp.getHeight();
		double[] meanProfile = new double[widthOut];
		double[] profile = new double[widthOut];
		double[] correlMean = new double[heightOut];
		
		for (int i = 0; i < widthOut; i++) {
			rotImp.setRoi(new Rectangle(i, 0, 1, heightOut));
			ImageStatistics is = rotImp.getStatistics(ImageStatistics.MEAN);
			meanProfile[i] = is.mean;
		}
		rotImp.deleteRoi();
		ImageProcessor rotProc = rotImp.getProcessor();
		for (int j = 0; j < heightOut; j++) {
			for (int i = 0; i < widthOut; i++) {
				profile[i] = rotProc.getPixelValue(i, j);
			}
			correlMean[j] = getCorrelation2(meanProfile, profile);
		}
		rotImp.changes = false;
		rotImp.close();
		
		double[] correlMeanAndSD = getMeanAndSD(correlMean);
		return correlMeanAndSD[0];
	}
	
	private double getCorrelation2(double[] arr1, double[] arr2) {
		double[] arr1MeanSD = getMeanAndSD(arr1);
		double[] arr2MeanSD = getMeanAndSD(arr2);
		double res = 0;
		for (int i = 0; i < arr1.length; i++) {
			res += (arr1[i]-arr1MeanSD[0]) * (arr2[i]-arr2MeanSD[0]);
		}
		res /= (arr1.length * arr1MeanSD[1] * arr2MeanSD[1]);
		return res;
	}
		
	private double[] getMeanAndSD(double[] arr) {
		double[] res = new double[2];
		res[0] = 0; res[1] = 0;
		int count = 0;
		for (int i = 0; i < arr.length; i++) {
			if (!Double.isNaN(arr[i])) {
				res[0] += arr[i];
				count += 1;
			}
		}
		res[0] /= count;
		count = 0;
		for (int i = 0; i < arr.length; i++) {
			if (!Double.isNaN(arr[i])) {
				res[1] += Math.pow(arr[i] - res[0], 2);
				count += 1;
			}
		}
		res[1] = Math.sqrt(res[1] / count);
		return res;
	}
	
	private double computeNormalizedVerticalVariance(ImagePlus imp, double maxVar) {
		int width = imp.getWidth();
		int height = imp.getHeight();
		varianceArray = new double[width];
		double maxRes = Double.NEGATIVE_INFINITY;
		for (int iWidth = 0; iWidth < width; iWidth++) {
			imp.setRoi(iWidth, 0, 1, height);
			ImageStatistics is = imp.getStatistics(ImageStatistics.STD_DEV);
			varianceArray[iWidth] = is.stdDev * is.stdDev;
			if (varianceArray[iWidth] > maxRes) {
				maxRes = varianceArray[iWidth];
			}
		}
		if (maxVar > 0)
			maxRes = maxVar;
		for (int iWidth = 0; iWidth < width; iWidth++) {
			varianceArray[iWidth] /= maxRes;
		}
		IJ.run(imp, "Select None", "");
		return maxRes;
	}
	
	private void setColorGroup(int groupId) {
		if (getRM().getSelectedIndex() > -1) {
			Roi[] r = getRM().getSelectedRoisAsArray();
			for (int iRoi = 0; iRoi < r.length; iRoi++) {
				r[iRoi].setStrokeColor(Colors.decode(ROI_COLORS[groupId]));
			}
		} else if (IJ.getImage().getRoi() != null) {
			IJ.getImage().getRoi().setStrokeColor(Colors.decode(ROI_COLORS[groupId]));
		}
	}
	
	private int getColorGroup(Roi r) {
		int res = 0;
		String roiColorStr = Colors.colorToString(r.getStrokeColor());
		for (int iColor = 0; iColor < ROI_COLORS.length; iColor++) {
			if (roiColorStr.matches(ROI_COLORS[iColor])) {
				res = iColor;
				break;
			}
		}
		return res;
	}
	
	public static double secondOrderPolynomialMax(double[] xArray, double[] yArray) {
		// Using two arrays of x, y coordinates, find the maximum of y, then compute
		// the parabolla that intersect the maximum and its two neighbours. Return
		// the x coordinate of this max.
		int res = 0;
		if (xArray.length != yArray.length) {
		    return res;
		}
		int maxIdx = 0;
		double yArrayMax = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < yArray.length; i++) {
			if (yArray[i] > yArrayMax) {
				yArrayMax = yArray[i];
				maxIdx = i;
			}
		}
		double x1, x2, x3, y1, y2, y3;
		if (maxIdx == 0) {
		    x1 = 2*xArray[maxIdx]-xArray[maxIdx+1]; // extrapol value for x (assuming non sparsity)
		    y1 = yArray[yArray.length-1]; // use opposite end value for y (assuming circularity)
		} else {
		    x1 = xArray[maxIdx-1];
		    y1 = yArray[maxIdx-1];
		}
		x2 = xArray[maxIdx];
		y2 = yArray[maxIdx];
		if (maxIdx == xArray.length-1) {
		    x3 = 2*xArray[maxIdx]-xArray[maxIdx-1]; // extrapol value for x (assuming non sparsity)
		    y3 = yArray[0]; // use opposite end value for y (assuming circularity)
		} else {
		    x3 = xArray[maxIdx+1];
		    y3 = yArray[maxIdx+1];
		}

		double a = ((y3-y1)/(x3-x1)-(y2-y1)/(x2-x1))/((x3*x3-x1*x1)/(x3-x1)-(x2*x2-x1*x1)/(x2-x1));
		double b = (y2-y1-a*(x2*x2-x1*x1))/(x2-x1);
		//double c = y1-a*x1*x1-b*x1;
		return -b/(2*a);
	}
	
	private ImagePlus pretreatXTImage(ImagePlus imp, String method, boolean inPlace) {
		ImagePlus res = new ImagePlus();
		if (method.matches("Demeaning"))
			res = imageDemeaning(imp, inPlace);
		else if (method.matches("SobelX2"))
			res = imageConvolve(imp,"Sobel2",false, inPlace);
		else if (method.matches("SobelXY2"))
			res = imageConvolve(imp,"Sobel2",true, inPlace);
		else if (method.matches("SobelX3"))
			res = imageConvolve(imp,"Sobel3",false, inPlace);
		else if (method.matches("SobelXY3"))
			res = imageConvolve(imp,"Sobel3",true, inPlace);
		else if (method.matches("SobelX5"))
			res = imageConvolve(imp,"Sobel5",false, inPlace);
		else if (method.matches("SobelXY5"))
			res = imageConvolve(imp,"Sobel5",true, inPlace);
		else if (method.matches("PrewittX3"))
			res = imageConvolve(imp,"Prewitt3",false, inPlace);
		else if (method.matches("PrewittXY3"))
			res = imageConvolve(imp,"Prewitt3",true, inPlace);
		return res;
	}
	
	private ImagePlus imageDemeaning(ImagePlus imp, boolean inPlace) {
		int width = imp.getWidth();
		int height = imp.getHeight();
		ImagePlus res;
		if (inPlace) {
			if (imp.getType() != ImagePlus.GRAY32)
				imp.setProcessor(imp.getProcessor().convertToFloat());
			res = imp;
		} else {
			res = imp.duplicate();
			res.setTitle("Demeaning");
			if (res.getType() != ImagePlus.GRAY32)
				res.setProcessor(res.getProcessor().convertToFloat());
		}
		for (int iHeight = 0; iHeight < height; iHeight++) {
			res.setRoi(0, iHeight, width, 1);
			ImageStatistics is = res.getStatistics(ImageStatistics.MEAN + ImageStatistics.STD_DEV);
			IJ.run(res, "Subtract...", "value="+is.mean);
			IJ.run(res, "Divide...", "value="+is.stdDev);
		}
		IJ.run(res, "Select None", "");
		IJ.resetMinAndMax(res);
		return res;
	}
	
	private ImagePlus imageConvolve(ImagePlus imp, String filterName, boolean bothDirections, boolean inPlace) {
		String dirStr = "X";
		if (bothDirections) {
			dirStr = "XY";
		}
		ImagePlus res;
		if (inPlace) {
			res = imp;
		} else {
			res = imp.duplicate();
			res.setTitle(filterName+dirStr);
		}
		String horiFilter = "[1]";
		String vertFilter = "[1]";
		if (filterName == "Sobel2") {
			horiFilter = "[1 -1 0\n1 -1 0\n0 0 0]";
			vertFilter = "[1 1 0\n-1 -1 0\n0 0 0]";
		}
		if (filterName == "Sobel3") {
			horiFilter = "[1 0 -1\n2 0 -2\n1 0 -1]";
			vertFilter = "[1 2 1\n0 0 0\n-1 -2 -1]";
		}
		if (filterName == "Sobel5") {
			horiFilter = "[5 4 0 -4 -5\n8 10 0 -10 -8\n10 20 0 -20 -10\n8 10 0 -10 -8\n5 4 0 -4 -5]";
			vertFilter = "[5 8 10 8 5\n4 10 20 10 4\n0 0 0 0 0\n-4 -10 -20 -10 -4\n-5 -8 -10 -8 -5]";
		}
		if (filterName == "Prewitt3") {
			horiFilter = "[1 0 -1\n1 0 -1\n1 0 -1]";
			vertFilter = "[1 1 1\n0 0 0\n-1 -1 -1]";
		}
		ImagePlus resVert = new ImagePlus();
		if (bothDirections) {
			resVert = imp.crop();
			IJ.run(resVert, "Convolve...", "text1="+vertFilter+" normalize");
			IJ.run(resVert, "Square", "");
		}
		IJ.run(res, "Convolve...", "text1="+horiFilter+" normalize");
		if (bothDirections) {
			IJ.run(res, "Square", "");
			res = ImageCalculator.run(res, resVert, "Add");
			IJ.run(res, "Square Root", "");
			resVert.close();
		}
		return res;
	}
	
	private void manualVerification() {
		if (emptyListAbort()) return;
		getVerif();
		nextVerification();
	}
	
	private void updateAngle() {
		int rtLine = getRTLineFromROI();
		if (rtLine == -1) return;
		Point[] p = verifImp.getRoi().getContainedPoints();
		double newAngle = (Math.atan2(p[p.length-1].getX()-p[0].getX(), p[p.length-1].getY()-p[0].getY())*180.0/Math.PI)%180;
		if (newAngle < 0) newAngle += 180;
		double newVelocity = rt.getValue(YPX_UM_COL_NAME,rtLine)*1000/(rt.getValue(XPX_MS_COL_NAME,rtLine)*Math.tan(newAngle*Math.PI/180.0));
		rt.setValue(VELOCITY_OK_COL_NAME, rtLine, newVelocity);
		rt.setValue(ANGLE_OK_COL_NAME, rtLine, newAngle);
		rt.setValue(VERIF_COL_NAME, rtLine, 1);
		rt.save(rtManualPath);
		rt.show(VELOCITY_MANUAL_FILENAME);
		updateStripes(verifImp, newAngle, getBoundaries());
	}
	
	private void discardAngle() {
		int rtLine = getRTLineFromROI();
		if (rtLine == -1) return;
		rt.setValue(VELOCITY_OK_COL_NAME, rtLine, Double.NaN);
		rt.setValue(ANGLE_OK_COL_NAME, rtLine, Double.NaN);
		rt.setValue(VERIF_COL_NAME, rtLine, -1);
		rt.save(rtManualPath);
		rt.show(VELOCITY_MANUAL_FILENAME);
		updateStripes(verifImp, Double.NaN, getBoundaries());
	}
	
	private void openStackVerif() {
		int rtLine = getRTLineFromROI();
		if (rtLine == -1) return;
		if (autoCloseCB.getState()) {
			if (oriImp != null) oriImp.close();
		}
		stacksListValue.select(rt.getStringValue(FILE_COL_NAME, rtLine));
		openStack(openVirtualCB.getState(), false);
		getRM().select((int)rt.getValue(ROI_COL_NAME, rtLine)-1);
	}
	
	private void openXTRawVerif() {
		int rtLine = getRTLineFromROI();
		if (rtLine == -1) return;
		if (autoCloseCB.getState()) {
			if (xtImp != null) xtImp.close();
		}
		setPaths(directory+rt.getStringValue(FILE_COL_NAME, rtLine));
		String xtFullPath = xtPath+"XT"+rt.getStringValue(ROI_COL_NAME, rtLine)+"_"+oriName+".tif";
		if (new File(xtFullPath).exists()) {
			xtImp = IJ.openImage(xtFullPath);
			xtImp.show();
		}
	}
	
	private void validateROI() {
		if (idxUID+1 >= rangeUID.length) {
			new MessageDialog(new Frame(), "Done", "Verification done.");
			return;
		}
		for (int rtLine = rangeUID[idxUID]; rtLine < rangeUID[idxUID+1]; rtLine++) {
			if (rt.getValue(VERIF_COL_NAME,rtLine) == 0) {
				if (rt.getValue(GROUP_COL_NAME, rtLine) == 3) {
					rt.setValue(VERIF_COL_NAME, rtLine, -1);
					rt.setValue(VELOCITY_OK_COL_NAME, rtLine, Double.NaN);
					rt.setValue(ANGLE_OK_COL_NAME, rtLine, Double.NaN);
				} else {
					rt.setValue(VERIF_COL_NAME, rtLine, 2);
					rt.setValue(VELOCITY_OK_COL_NAME, rtLine, rt.getValue(VELOCITY_COL_NAME, rtLine));
					rt.setValue(ANGLE_OK_COL_NAME, rtLine, rt.getValue(ANGLE_COL_NAME, rtLine));
				}
			}
		}
		rt.save(rtManualPath);
		rt.show(VELOCITY_MANUAL_FILENAME);
		nextVerification();
	}
	
	private void unvalidateROI() {
		if (idxUID+1 >= rangeUID.length) {
			new MessageDialog(new Frame(), "Done", "Verification done.");
			return;
		}
		for (int rtLine = rangeUID[idxUID]; rtLine < rangeUID[idxUID+1]; rtLine++) {
			rt.setValue(VERIF_COL_NAME, rtLine, -2);
			rt.setValue(VELOCITY_OK_COL_NAME, rtLine, Double.NaN);
			rt.setValue(ANGLE_OK_COL_NAME, rtLine, Double.NaN);
		}
		rt.save(rtManualPath);
		rt.show(VELOCITY_MANUAL_FILENAME);
		nextVerification();
	}
	
	private int getRTLineFromROI() {
		int[] boundaries = getBoundaries();
		if (boundaries == null) {
			return -1;
		} else {
			int offsetX = boundaries[0]/boundaries[2];
			int offsetY = boundaries[1]/boundaries[3];
			int offsetIndex = offsetX+xyUID[0]*offsetY;
			int indexLimit;
			if (idxUID+1 >= rangeUID.length) {
				indexLimit = rangeUID.length;
			} else {
				indexLimit = rangeUID[idxUID+1];
			}
			if (offsetIndex+rangeUID[idxUID] >= indexLimit) {
				return -1;
			} else {
				return offsetIndex+rangeUID[idxUID];
			}
		}
	}
	
	private int[] getBoundaries() {
		if (verifImp == null) {
			new MessageDialog(new Frame(), "Image missing", "XT Verification image is not open.");
			return null;
		}
		Roi r = verifImp.getRoi();
		if (r == null || r.getType() != 5) {
			new MessageDialog(new Frame(), "Line ROI missing", "A line ROI is needed to update values.");
			return null;
		}
		int[] res = new int[4];
		Rectangle rect = r.getBounds();
		double xPos = rect.x + rect.width/2;
		double yPos = rect.y + rect.height/2;
		res[0] = (int)Math.floor(xPos / whUID[0]) * whUID[0];
		res[1] = (int)Math.floor(yPos / whUID[1]) * whUID[1];
		res[2] = whUID[0] - BORDER_SIZE;
		res[3] = whUID[1] - BORDER_SIZE;
		return res;
	}
	
	private void nextVerification() {
		rtPath = directory+VELOCITY_FILENAME;
		rtManualPath = directory+VELOCITY_MANUAL_FILENAME;
		if (rtPath == null || rtManualPath == null) return;
		File rtFile = new File(rtManualPath);
		boolean openManual = true;
		if (!rtFile.exists()) {
			rtFile = new File(rtPath);
			if (!rtFile.exists()) {
				new MessageDialog(new Frame(), "Missing table", "No results table to verify.");
				return;
			}
			openManual = false;
		}
		try {
			if (openManual) {
				rt = ResultsTable.open(rtManualPath);
			} else {
				rt = ResultsTable.open(rtPath);
				rt.setValue(VERIF_COL_NAME, 0, 0);
			}
			rt.setPrecision(TABLE_PRECISION);
			rt.showRowNumbers(false);
		} catch (IOException e) {
			if (openManual) {
				new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtManualPath);
			} else {
				new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtPath);
			}
			return;
		}
		rt.show(VELOCITY_MANUAL_FILENAME);
		rangeUID = new int[rt.size()+1];
		rangeUID[0] = 0;
		idxUID = 1;
		for (int iRT = 1; iRT < rt.size(); iRT++) {
			String prevID = "" + rt.getValue(ROI_COL_NAME, iRT-1); // + rt.getStringValue(PREPROCESS_COL_NAME, iRT-1);
			String currID = "" + rt.getValue(ROI_COL_NAME, iRT); // + rt.getStringValue(PREPROCESS_COL_NAME, iRT);
			if (currID.compareTo(prevID) != 0) {
				rangeUID[idxUID++] = iRT;
			}
		}
		rangeUID[idxUID++] = rt.size();
		rangeUID = Arrays.copyOf(rangeUID, idxUID);
		boolean hasUnmodifiedRows = false;
		for (idxUID = 0; idxUID < rangeUID.length-1; idxUID++) {
			for (int rtLine = rangeUID[idxUID]; rtLine < rangeUID[idxUID+1]; rtLine++) {
				if (rt.getValue(VERIF_COL_NAME,rtLine) == 0) {
					hasUnmodifiedRows = true;
					break;
				}
			}
			if (hasUnmodifiedRows) break;
		}
		int idxImage = 0;
		if (verifImp != null) verifImp.close();
		if (autoCloseCB.getState()) {
			if (oriImp != null) oriImp.close();
			if (xtImp != null) xtImp.close();
		}
		if (idxUID == rangeUID.length-1) {
			new MessageDialog(new Frame(), "Done", "Verification done.");
			return;
		}
		verifImp = new ImagePlus();
		for (int rtLine = rangeUID[idxUID]; rtLine < rangeUID[idxUID+1]; rtLine++) {
			String name = rt.getStringValue(FILE_COL_NAME, rtLine);
			File f = new File(directory+name);
			name = f.getName();
			int lastIdx = name.lastIndexOf(".");
			name = name.substring(0, lastIdx)+".tif";
			String xtFullPath = f.getParent()+File.separator+"XT"+File.separator+"XT"+(int)rt.getValue(ROI_COL_NAME, rtLine)+"_"+name;
			ImagePlus tmpImp = IJ.openImage(xtFullPath);
			if (tmpImp.getType() != ImagePlus.GRAY32) {
				tmpImp.setProcessor(tmpImp.getProcessor().convertToFloat());
			}
			//tmpImp = pretreatXTImage(tmpImp, rt.getStringValue(PREPROCESS_COL_NAME, rtLine), true);
			tmpImp = pretreatXTImage(tmpImp, "SobelX3", true);
			
			double angle = rt.getValue(ANGLE_COL_NAME, rtLine);
			if (rt.getValue(VERIF_COL_NAME, rtLine) != 0) {
				angle = rt.getValue(ANGLE_OK_COL_NAME, rtLine);
			}
			double[] angles;
			if (Double.isNaN(angle)) {
				angles = new double[2];
				angles[0] = -45;
				angles[1] = 45;
			} else {
				angles = new double[1];
				angles[0] = angle;
			}
			int[] boundaries = {0,0,tmpImp.getWidth(),tmpImp.getHeight()};
			plotWithStripes(tmpImp, angles, 0, boundaries);
			tmpImp = tmpImp.flatten();
			if (idxImage == 0) {
				verifImp = tmpImp;
			} else {
				verifImp = Concatenator.run(verifImp, tmpImp);
			}
			idxImage++;
		}
		whUID[0] = verifImp.getWidth()+BORDER_SIZE;
		whUID[1] = verifImp.getHeight()+BORDER_SIZE;
		xyUID = getGrid(idxImage, whUID[0], whUID[1]);
		IJ.run(verifImp, "Make Montage...", "columns="+xyUID[0]+" rows="+xyUID[1]+" scale=1 border="+BORDER_SIZE+" label");
		verifImp.close();
		verifImp = IJ.getImage();
		verifImp.setTitle(ROI_COL_NAME+rt.getStringValue(ROI_COL_NAME,rangeUID[idxUID]));
		// verifImp.setTitle(ROI_COL_NAME+rt.getStringValue(ROI_COL_NAME,rangeUID[idxUID])+"-"+rt.getStringValue(PREPROCESS_COL_NAME,rangeUID[idxUID]));
		IJ.setTool(4);
	}
	
	private void updateStripes(ImagePlus verifImp, double newAngle, int[] boundaries) {
		double[] angles;
		if (Double.isNaN(newAngle)) {
			angles = new double[2];
			angles[0] = -45;
			angles[1] = 45;
		} else {
			angles = new double[1];
			angles[0] = newAngle;
		}
		removeStripes(verifImp, boundaries);
		plotWithStripes(verifImp, angles, 1, boundaries);
	}
	
	private void plotWithStripes(ImagePlus imp, double[] angles, int colorIdx, int[] boundaries) {
		double stripeDist = 100; // target distance between stripes
		int width = boundaries[2]; //imp.getWidth();
		int height = boundaries[3]; //imp.getHeight();
		if (colorIdx == 0) {
			IJ.resetMinAndMax(imp);
			IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		}
		int[] colors = {1730252,15093018,15119130,8395417,5092122,1753074,10033690}; // Matlab plot colors
		int color = colors[colorIdx % colors.length];
		for (int iAngle = 0; iAngle < angles.length; iAngle++) {
			double theta = angles[iAngle];
			double d = getHorizontalMaxDistance(width, height, theta);
			double offset = height/2.0*Math.tan(theta*Math.PI/180);
			
			int nStripe = (int)(Math.round(width*Math.abs(Math.cos(theta*Math.PI/180))/2/stripeDist)*2+1);
			double step = 2.0/nStripe;
			double start = 1.0-step/2.0;
			for (double j = -start*d+width/2.0; j <= start*d+width/2.0; j+=step*d) {
				double xStart = j-offset;
				double yStart = 0;
				double xStop = j+offset;
				double yStop = height;
				if (xStart < 0) {
					yStart = - xStart / Math.tan(theta*Math.PI/180);
					xStart = 0;
				} else if (xStart > width) {
					yStart = (width - xStart) / Math.tan(theta*Math.PI/180);
					xStart = width;
				}
				if (xStop > width) {
					yStop = height - (width - xStop) / Math.tan((180-theta)*Math.PI/180);
					xStop = width;
				} else if (xStop < 0) {
					yStop = height + xStop / Math.tan((180-theta)*Math.PI/180);
					xStop = 0;
				}
				imp.setRoi(new Line(xStart+boundaries[0], yStart+boundaries[1], xStop+boundaries[0], yStop+boundaries[1]));
				Roi r = imp.getRoi();
				r.setStrokeColor(new Color(color));
				r.setStrokeWidth(2);
				IJ.run(imp, "Add Selection...", "");
			}
		}
		IJ.run(imp, "Select None", "");
	}
	
	private void removeStripes(ImagePlus imp, int[] boundaries) {
		Overlay o = imp.getOverlay();
		if (o != null) {
			for (int iRoi = o.size()-1; iRoi >= 0; iRoi--) {
				Rectangle r = o.get(iRoi).getBounds();
				double centerX = (double)r.x + (double)r.width / 2.0;
				double centerY = (double)r.y + (double)r.height / 2.0;
				if (centerX >= boundaries[0] && centerX <= boundaries[0] + boundaries[2] && centerY >= boundaries[1] && centerY <= boundaries[1] + boundaries[3]) {
					o.remove(iRoi);
				}
			}
		}
	}
	
	private double getHorizontalMaxDistance(int width, int height, double theta) {
		// Finds the highest horizontal distance between a line theta degree from vertical and that
		// crosses a width * height rectangle and the center of that rectangle
		double dist;
	    if (theta == 0 || theta == 180) {
	        dist = width/2;
	    } else {
	        // Equation of the line theta degrees from the vertical :
	        // a*x + b*y + c = 0;
	        // Distance between the line and the point (x0, y0) :
	        // dist = |a*x0+b*y0+c| / d, with d = sqrt(a*a+b*b)
	        double a = 1/Math.tan(theta*Math.PI/180);
	        double b = -1;
	        double c = (height - a * width)/2;
	        double d = Math.sqrt(a*a+b*b);
	        double d1 = Math.abs(c)/d;
	        double d2 = Math.abs(a*width+c)/d;
	        double d3 = Math.abs(b*height+c)/d;
	        double d4 = Math.abs(a*width+b*height+c)/d;
	        dist = Math.abs(Math.max(Math.max(d1,d2),Math.max(d3,d4))/Math.cos(theta*Math.PI/180));
	    }
	    return dist;
	}
	
	private int[] getGrid(int nSlices, int width, int height) {
		int[] res = new int[2];
		res[0] = nSlices; res[1] = 1;
		double targetRatio = screenWidth/screenHeight;
		double prevRatio = nSlices * width / height;
		for (int i = 2; i <= nSlices; i++) {
			double y = i;
			double x = Math.ceil(nSlices/y);
			double newRatio = (x*width) / (y*height);
			if (Math.abs(targetRatio-prevRatio) > Math.abs(targetRatio-newRatio)) {
				prevRatio = newRatio;
				res[1] = i;
				res[0] = (int)Math.ceil(nSlices/y);
				if ((res[1]-1)*res[0] >= nSlices) {
					res[1] -= 1;
				}
			}
		}
		return res;
	}
	
	private void evaluateVerification() {
		if (emptyListAbort()) return;
		if (rtManualPath == null) return;
		File rtFile = new File(rtManualPath);
		if (!rtFile.exists()) {
			new MessageDialog(new Frame(), "Missing table", "No results table to verify.");
			return;
		}
		try {
			rt = ResultsTable.open(rtManualPath);
			rt.setPrecision(TABLE_PRECISION);
			rt.showRowNumbers(false);
		} catch (IOException e) {
			new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtManualPath);
			return;
		}
		rt.show(VELOCITY_MANUAL_FILENAME);
		
		GenericDialog gd = new GenericDialog("Verification metrics");
		//gd.addMessage("Global parameters");
		gd.addChoice("Metric 1: ", VERIF_METRICS, Prefs.get("RBCv.verifMetric1", VERIF_METRICS[0]));
		gd.addChoice("Metric 2: ", VERIF_METRICS, Prefs.get("RBCv.verifMetric2", VERIF_METRICS[2]));
		gd.showDialog();
		if (gd.wasCanceled()) return;
		String verifMetric1 = gd.getNextChoice();
		Prefs.set("RBCv.verifMetric1", verifMetric1);
		String verifMetric2 = gd.getNextChoice();
		Prefs.set("RBCv.verifMetric2", verifMetric2);
		double[] metric1All = rt.getColumn(verifMetric1);
		double[] metric2All = rt.getColumn(verifMetric2);
		boolean m1above = false, m2above = false;
		for (int iVerif = 0; iVerif < VERIF_METRICS.length; iVerif++) {
			if (verifMetric1.matches(VERIF_METRICS[iVerif])) m1above = VERIF_METRICS_ABOVE[iVerif];
			if (verifMetric2.matches(VERIF_METRICS[iVerif])) m2above = VERIF_METRICS_ABOVE[iVerif];
		}
		
		double[] verified = rt.getColumn(VERIF_COL_NAME);
		int[] verifGood = new int[verified.length];
		int[] verifBad = new int[verified.length];
		int[] verifNothing = new int[verified.length];
		int nGood = 0, nBad = 0, nNothing = 0;
		for (int rtLine = 0; rtLine < verified.length; rtLine++) {
			if (verified[rtLine] == 2) {
				verifGood[nGood++] = rtLine;
			} else if (verified[rtLine] == 1) {
				verifBad[nBad++] = rtLine;
			} else if (verified[rtLine] == -1) {
				verifNothing[nNothing++] = rtLine;
			}
		}
		double[] metric1Good = new double[nGood];
		double[] metric2Good = new double[nGood];
		double[] metric1Bad = new double[nBad];
		double[] metric2Bad = new double[nBad];
		double[] metric1Nothing = new double[nNothing];
		double[] metric2Nothing = new double[nNothing];
		double[] metric1Negative = new double[nBad+nNothing];
		double[] metric2Negative = new double[nBad+nNothing];
		for (int iGood = 0; iGood < nGood; iGood++) {
			metric1Good[iGood] = metric1All[verifGood[iGood]];
			metric2Good[iGood] = metric2All[verifGood[iGood]];
		}
		for (int iBad = 0; iBad < nBad; iBad++) {
			metric1Bad[iBad] = metric1All[verifBad[iBad]];
			metric2Bad[iBad] = metric2All[verifBad[iBad]];
			metric1Negative[iBad] = metric1All[verifBad[iBad]];
			metric2Negative[iBad] = metric2All[verifBad[iBad]];
		}
		for (int iNothing = 0; iNothing < nNothing; iNothing++) {
			metric1Nothing[iNothing] = metric1All[verifNothing[iNothing]];
			metric2Nothing[iNothing] = metric2All[verifNothing[iNothing]];
			metric1Negative[iNothing+nBad] = metric1All[verifNothing[iNothing]];
			metric2Negative[iNothing+nBad] = metric2All[verifNothing[iNothing]];
		}

		double min1 = Double.POSITIVE_INFINITY;
		double min2 = Double.POSITIVE_INFINITY;
		double max1 = Double.NEGATIVE_INFINITY;
		double max2 = Double.NEGATIVE_INFINITY;
		for (int iMeasure = 0; iMeasure < metric1All.length; iMeasure++) {
			if (metric1All[iMeasure] == metric1All[iMeasure]) { // to avoid NaNs
				min1 = Math.min(min1, metric1All[iMeasure]);
				max1 = Math.max(max1, metric1All[iMeasure]);
			}
			if (metric2All[iMeasure] == metric2All[iMeasure]) { // to avoid NaNs
				min2 = Math.min(min2, metric2All[iMeasure]);
				max2 = Math.max(max2, metric2All[iMeasure]);
			}
		}
		
		double step1 = Math.pow(10,Math.round(Math.log10((max1-min1)/100)));
		double step2 = Math.pow(10,Math.round(Math.log10((max2-min2)/100)));
		min1 = step1 * Math.floor(min1 / step1);
		max1 = step1 * Math.ceil(max1 / step1);
		min2 = step2 * Math.floor(min2 / step2);
		max2 = step2 * Math.ceil(max2 / step2);
		int n1 = (int)((max1-min1)/step1)+1;
		int n2 = (int)((max2-min2)/step2)+1;
		float[][] fnArr = new float[n1][n2];
		float[][] tnArr = new float[n1][n2];
		for (int iGood = 0; iGood < nGood; iGood++) {
			fnArr[(int)((metric1Good[iGood]-min1)/step1)][(int)((metric2Good[iGood]-min2)/step2)] += 1;
		}
		for (int iNeg = 0; iNeg < metric1Negative.length; iNeg++) {
			tnArr[(int)((metric1Negative[iNeg]-min1)/step1)][(int)((metric2Negative[iNeg]-min2)/step2)] += 1;
		}
		if (m1above) {
			for (int i1 = 1; i1 < n1; i1++) {
				for (int i2 = 0; i2 < n2; i2++) {
					fnArr[i1][i2] += fnArr[i1-1][i2];
					tnArr[i1][i2] += tnArr[i1-1][i2];
				}
			}
		} else {
			for (int i1 = n1-2; i1 >= 0; i1--) {
				for (int i2 = 0; i2 < n2; i2++) {
					fnArr[i1][i2] += fnArr[i1+1][i2];
					tnArr[i1][i2] += tnArr[i1+1][i2];
				}
			}
		}
		if (m2above) {
			for (int i1 = 0; i1 < n1; i1++) {
				for (int i2 = 1; i2 < n2; i2++) {
					fnArr[i1][i2] += fnArr[i1][i2-1];
					tnArr[i1][i2] += tnArr[i1][i2-1];
				}
			}
		} else {
			for (int i1 = 0; i1 < n1; i1++) {
				for (int i2 = n2-2; i2 >= 0; i2--) {
					fnArr[i1][i2] += fnArr[i1][i2+1];
					tnArr[i1][i2] += tnArr[i1][i2+1];
				}
			}
		}
		
		ImagePlus fnImp = IJ.createImage("fn", "32-bit black", n1, n2, 1);
		fnImp.getProcessor().setFloatArray(fnArr);
		ImagePlus tnImp = IJ.createImage("tn", "32-bit black", n1, n2, 1);
		tnImp.getProcessor().setFloatArray(tnArr);
		ImagePlus tpImp = fnImp.crop();
		ImageStatistics is = tpImp.getStatistics(ImageStatistics.MIN_MAX);
		IJ.run(tpImp, "Subtract...", "value="+is.max);
		IJ.run(tpImp, "Multiply...", "value=-1");
		ImagePlus fpImp = tnImp.crop();
		is = fpImp.getStatistics(ImageStatistics.MIN_MAX);
		IJ.run(fpImp, "Subtract...", "value="+is.max);
		IJ.run(fpImp, "Multiply...", "value=-1");
		// MCC = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
		ImagePlus tptnImp = ImageCalculator.run(tpImp, tnImp, "Multiply create 32-bit");
		ImagePlus fpfnImp = ImageCalculator.run(fpImp, fnImp, "Multiply create 32-bit");
		ImageCalculator.run(tptnImp, fpfnImp, "Subtract");
		ImagePlus tpfpImp = ImageCalculator.run(tpImp, fpImp, "Add create 32-bit");
		ImagePlus tpfnImp = ImageCalculator.run(tpImp, fnImp, "Add create 32-bit");
		ImagePlus tnfpImp = ImageCalculator.run(tnImp, fpImp, "Add create 32-bit");
		ImagePlus tnfnImp = ImageCalculator.run(tnImp, fnImp, "Add create 32-bit");
		ImageCalculator.run(tpfpImp, tpfnImp, "Multiply");
		ImageCalculator.run(tpfpImp, tnfpImp, "Multiply");
		ImageCalculator.run(tpfpImp, tnfnImp, "Multiply");
		IJ.run(tpfpImp, "Square Root", "");
		ImagePlus mccImp = ImageCalculator.run(tptnImp, tpfpImp, "Divide create 32-bit");
		mccImp.setTitle("Matthews correlation coefficient");
		IJ.run(mccImp, "Properties...", "pixel_width="+step1+" pixel_height="+(-step2)+" origin="+(-min1/step1)+","+(-min2/step2));
		mccImp.getCalibration().setXUnit(verifMetric1);
		mccImp.getCalibration().setYUnit(verifMetric2);
		IJ.run(mccImp, "Flip Vertically", "");
		mccImp.show();
		IJ.run(mccImp, "Find Maxima...", "prominence=1 output=[Point Selection]");
		Roi r = mccImp.getRoi();
		Point[] p = new Point[0];
		if (r != null) {
			p = r.getContainedPoints();
		}
		double[] xCoord = new double[p.length];
		double[] yCoord = new double[p.length];
		double[] xCoordMeasure = new double[p.length];
		double[] yCoordMeasure = new double[p.length];
		for (int iPoint = 0; iPoint < p.length; iPoint++) {
			xCoord[iPoint] = p[iPoint].getX();
			yCoord[iPoint] = p[iPoint].getY();
			xCoordMeasure[iPoint] = xCoord[iPoint] * step1 + min1;
			yCoordMeasure[iPoint] = step2 * Math.round((max2 - yCoord[iPoint] * step2)/step2);
			//IJ.log("Max "+(iPoint+1)+" image coords ["+xCoord[iPoint]+","+yCoord[iPoint]+"], "+verifMetric1+": "+
			//			xCoordMeasure[iPoint]+", "+verifMetric2+": "+yCoordMeasure[iPoint]);
		}
		
		tptnImp.close();
		fpfnImp.close();
		tpfpImp.close();
		tpfnImp.close();
		tnfpImp.close();
		tnfnImp.close();
		boolean closeStuff = true;
		if (closeStuff) {
			fnImp.close();
			tnImp.close();
			tpImp.close();
			fpImp.close();
		} else {
			fnImp.show();
			tnImp.show();
			tpImp.setTitle("tp");
			tpImp.show();
			fpImp.setTitle("fp");
			fpImp.show();
		}
		
		Plot verifPlot = new Plot("Verification",verifMetric1,verifMetric2);
		verifPlot.setColor("green");
		verifPlot.add("cross", metric1Good, metric2Good);
		verifPlot.setColor("blue");
		verifPlot.add("cross", metric1Bad, metric2Bad);
		verifPlot.setColor("red");
		verifPlot.add("cross", metric1Nothing, metric2Nothing);
		String legend = "Good measure\tBad measure\tNothing to measure";
		verifPlot.setColor("black");
		for (int iPoint = 0; iPoint < p.length; iPoint++) {
			if (!Double.isNaN(mccImp.getProcessor().getPixelValue((int)p[iPoint].getX(), (int)p[iPoint].getY()))) {
				verifPlot.add("line", new double[]{m1above?min1:max1,xCoordMeasure[iPoint],xCoordMeasure[iPoint]},
						new double[]{yCoordMeasure[iPoint],yCoordMeasure[iPoint],m2above?min2:max2});
				legend = legend+"\t["+xCoordMeasure[iPoint]+", "+yCoordMeasure[iPoint]+"]: "+mccImp.getProcessor().getPixelValue((int)xCoord[iPoint],(int)yCoord[iPoint]);
			}
		}
		verifPlot.setLegend(legend, Plot.AUTO_POSITION);
		verifPlot.show();
		verifPlot.setLimits(min1, max1, min2, max2);
		Zoom.set(mccImp, (double)verifPlot.getSize().height / mccImp.getHeight());
		
	}
	
	private String[] getExistingTableNames(String str) {
		String[] rtPathsAll;
		String[] velocityFileNamesAll;
		if (str.matches("Raw")) {
			rtPathsAll = new String[] {rtPath,rtSlidingPath};
			velocityFileNamesAll = new String[] {VELOCITY_FILENAME, VELOCITY_SLIDING_FILENAME};
		} else {
			rtPathsAll = new String[] {rtPath,rtManualPath,rtAutoPath,rtSlidingPath,rtSlidingAutoPath};
			velocityFileNamesAll = new String[] {VELOCITY_FILENAME, VELOCITY_MANUAL_FILENAME, VELOCITY_AUTO_FILENAME,
				VELOCITY_SLIDING_FILENAME, VELOCITY_SLIDING_AUTO_FILENAME};
		}
		boolean[] pathExists = new boolean[5];
		int nPathExists = 0;
		for (int i = 0; i < rtPathsAll.length; i++) {
			File rtFile = new File(rtPathsAll[i]);
			pathExists[i] = rtFile.exists();
			if (pathExists[i]) nPathExists++;
		}
		String[] velocityFileNames = new String[nPathExists];
		if (nPathExists > 0) {
			int pathIndex = 0;
			for (int i = 0; i < rtPathsAll.length; i++) {
				if (pathExists[i]) {
					velocityFileNames[pathIndex++] = velocityFileNamesAll[i];
				}
			}
		}
		return velocityFileNames;
	}
	
	private void automaticVerification() {
		if (emptyListAbort()) return;
		if (rtPath == null) return;
		
		String[] velocityFileNames = getExistingTableNames("Raw");
		if (velocityFileNames.length == 0) {
			new MessageDialog(new Frame(), "No result table", "No result table found. Aborting.");
			return;
		}
		
		GenericDialog gd = new GenericDialog("Verification metrics");
		//gd.addMessage("Global parameters");
		gd.addChoice("Choose table:", velocityFileNames, velocityFileNames[0]);
		gd.addChoice("Metric 1:", VERIF_METRICS, Prefs.get("RBCv.verifMetric1", VERIF_METRICS[0]));
		gd.addNumericField("Threshold Metric 1:", Prefs.get("RBCv.verifThreshold1", 10));
		gd.addChoice("Metric 2:", VERIF_METRICS, Prefs.get("RBCv.verifMetric2", VERIF_METRICS[2]));
		gd.addNumericField("Threshold Metric 2:", Prefs.get("RBCv.verifThreshold2", 0.58));
		gd.showDialog();
		if (gd.wasCanceled()) return;
		String velocityFileName = gd.getNextChoice();
		String verifMetric1 = gd.getNextChoice();
		Prefs.set("RBCv.verifMetric1", verifMetric1);
		double verifThreshold1 = gd.getNextNumber();
		Prefs.set("RBCv.verifThreshold1", verifThreshold1);
		String verifMetric2 = gd.getNextChoice();
		Prefs.set("RBCv.verifMetric2", verifMetric2);
		double verifThreshold2 = gd.getNextNumber();
		Prefs.set("RBCv.verifThreshold2", verifThreshold2);
		boolean m1above = false, m2above = false;
		for (int iVerif = 0; iVerif < VERIF_METRICS.length; iVerif++) {
			if (verifMetric1.matches(VERIF_METRICS[iVerif])) m1above = VERIF_METRICS_ABOVE[iVerif];
			if (verifMetric2.matches(VERIF_METRICS[iVerif])) m2above = VERIF_METRICS_ABOVE[iVerif];
		}
		
		boolean isManual = velocityFileName.matches(VELOCITY_FILENAME);
		String rtPathChoice = isManual?rtPath:rtSlidingPath;
		String velocityFileNameResult = isManual?VELOCITY_AUTO_FILENAME:VELOCITY_SLIDING_AUTO_FILENAME;
		String rtPathResult = isManual?rtAutoPath:rtSlidingAutoPath;
		
		try {
			rt = ResultsTable.open(isManual?rtPath:rtSlidingPath);
			rt.setPrecision(TABLE_PRECISION);
			rt.showRowNumbers(false);
		} catch (IOException e) {
			new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtPathChoice);
			return;
		}
		
		for (int rtLine = 0; rtLine < rt.size(); rtLine++) {
			if ((m1above ^ rt.getValue(verifMetric1, rtLine) < verifThreshold1) || (m2above ^ rt.getValue(verifMetric2,  rtLine) < verifThreshold2)) {
				rt.setValue(VERIF_COL_NAME, rtLine, 3);
				rt.setValue(VELOCITY_OK_COL_NAME, rtLine, rt.getValue(VELOCITY_COL_NAME, rtLine));
				rt.setValue(ANGLE_OK_COL_NAME, rtLine, rt.getValue(ANGLE_COL_NAME, rtLine));
			} else {
				rt.setValue(VERIF_COL_NAME, rtLine, -3);
				rt.setValue(VELOCITY_OK_COL_NAME, rtLine, Double.NaN);
				rt.setValue(ANGLE_OK_COL_NAME, rtLine, Double.NaN);
			}
		}
		rt.show(velocityFileNameResult);
		rt.save(rtPathResult);
	}
	
	private void velocityLimit() {
		if (emptyListAbort()) return;
		if (rtPath == null) return;
		
		String[] velocityFileNames = getExistingTableNames("All");
		if (velocityFileNames.length == 0) {
			new MessageDialog(new Frame(), "No result table", "No result table found. Aborting.");
			return;
		}
		
		GenericDialog gd = new GenericDialog("Velocity limit");
		gd.addChoice("Choose table:", velocityFileNames, velocityFileNames[0]);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		String velocityFileName = gd.getNextChoice();
		
		String rtPathChoice;
		if (velocityFileName.matches(VELOCITY_FILENAME)) {
			rtPathChoice = rtPath;
		} else if (velocityFileName.matches(VELOCITY_MANUAL_FILENAME)) {
			rtPathChoice = rtManualPath;
		} else if (velocityFileName.matches(VELOCITY_AUTO_FILENAME)) {
			rtPathChoice = rtAutoPath;
		} else if (velocityFileName.matches(VELOCITY_SLIDING_FILENAME)) {
			rtPathChoice = rtSlidingPath;
		} else { // if (velocityFileName.matches(VELOCITY_SLIDING_AUTO_FILENAME)) {
			rtPathChoice = rtSlidingAutoPath;
		}
		try {
			rt = ResultsTable.open(rtPathChoice);
			rt.setPrecision(TABLE_PRECISION);
			rt.showRowNumbers(false);
		} catch (IOException e) {
			new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtPathChoice);
			return;
		}
		
		boolean verified = (rt.getColumnIndex(VERIF_COL_NAME) != ResultsTable.COLUMN_NOT_FOUND);
		if (rt.getColumnIndex(VELOCITY_OK_COL_NAME) == ResultsTable.COLUMN_NOT_FOUND) {
			rt.setColumn(VELOCITY_OK_COL_NAME, rt.getColumnAsVariables(VELOCITY_COL_NAME));
		}
		if (rt.getColumnIndex(ANGLE_OK_COL_NAME) == ResultsTable.COLUMN_NOT_FOUND) {
			rt.setColumn(ANGLE_OK_COL_NAME, rt.getColumnAsVariables(ANGLE_COL_NAME));
		}
		
		for (int rtLine = 0; rtLine < rt.size(); rtLine++) {
			if (rtLine == 0 || !rt.getStringValue(FILE_COL_NAME, rtLine).matches(rt.getStringValue(FILE_COL_NAME, rtLine-1))
					|| !rt.getStringValue(ROI_COL_NAME, rtLine).matches(rt.getStringValue(ROI_COL_NAME, rtLine-1))) {
				setPaths(directory+rt.getStringValue(FILE_COL_NAME, rtLine));
				openImageROIs(true, false);
			}
			if (verified && rt.getValue(VERIF_COL_NAME,rtLine) < 1) continue;
			double roiLength = rt.getValue(LENGTH_COL_NAME, rtLine);
			double angleTan = Math.tan(Math.PI/180.0*rt.getValue(ANGLE_OK_COL_NAME, rtLine));
			
			if (Math.abs(angleTan) * roiLength < 2) {
				rt.setValue(VERIF_COL_NAME, rtLine, -4);
				rt.setValue(VELOCITY_OK_COL_NAME, rtLine, Double.NaN);
				rt.setValue(ANGLE_OK_COL_NAME, rtLine, Double.NaN);
			}
		}
		rt.save(rtPathChoice);
		rt.show(velocityFileName);
	}
	
	private void rollingShutterAdjust() {
		if (emptyListAbort()) return;
		if (rtPath == null) return;
		
		String[] velocityFileNames = getExistingTableNames("All");
		if (velocityFileNames.length == 0) {
			new MessageDialog(new Frame(), "No result table", "No result table found. Aborting.");
			return;
		}
		
		GenericDialog gd = new GenericDialog("Rolling shutter adjustment");
		gd.addChoice("Choose table:", velocityFileNames, velocityFileNames[0]);
		gd.addMessage("Camera parameters");
		gd.addNumericField("Line time", Prefs.get("RBCv.camLineTime", 10.0), 5, 8, "us");
		gd.addNumericField("Y offset", Prefs.get("RBCv.camYOffset", 0), 0, 8, "");
		gd.addNumericField("Rolling shutter start line", Prefs.get("RBCv.rollingShutterStartLine", 0), 0, 8, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		String velocityFileName = gd.getNextChoice();
		double camLineTime = gd.getNextNumber();
		Prefs.set("RBCv.camLineTime", camLineTime);
		int camYOffset = (int)(gd.getNextNumber());
		Prefs.set("RBCv.camYOffset", camYOffset);
		int rollingShutterStartLine = (int)(gd.getNextNumber());
		Prefs.set("RBCv.rollingShutterStartLine", rollingShutterStartLine);
		
		String rtPathChoice;
		if (velocityFileName.matches(VELOCITY_FILENAME)) {
			rtPathChoice = rtPath;
		} else if (velocityFileName.matches(VELOCITY_MANUAL_FILENAME)) {
			rtPathChoice = rtManualPath;
		} else if (velocityFileName.matches(VELOCITY_AUTO_FILENAME)) {
			rtPathChoice = rtAutoPath;
		} else if (velocityFileName.matches(VELOCITY_SLIDING_FILENAME)) {
			rtPathChoice = rtSlidingPath;
		} else { // if (velocityFileName.matches(VELOCITY_SLIDING_AUTO_FILENAME)) {
			rtPathChoice = rtSlidingAutoPath;
		}
		try {
			rt = ResultsTable.open(rtPathChoice);
			rt.setPrecision(TABLE_PRECISION);
			rt.showRowNumbers(false);
		} catch (IOException e) {
			new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtPathChoice);
			return;
		}
		
		boolean verified = (rt.getColumnIndex(VERIF_COL_NAME) != ResultsTable.COLUMN_NOT_FOUND);
		if (rt.getColumnIndex(VELOCITY_OK_COL_NAME) == ResultsTable.COLUMN_NOT_FOUND) {
			rt.setColumn(VELOCITY_OK_COL_NAME, rt.getColumnAsVariables(VELOCITY_COL_NAME));
		}
		if (rt.getColumnIndex(ANGLE_OK_COL_NAME) == ResultsTable.COLUMN_NOT_FOUND) {
			rt.setColumn(ANGLE_OK_COL_NAME, rt.getColumnAsVariables(ANGLE_COL_NAME));
		}
		
		double camLineTimeSec = camLineTime / 1000000;
		
		for (int rtLine = 0; rtLine < rt.size(); rtLine++) {
			if (rtLine == 0 || !rt.getStringValue(FILE_COL_NAME, rtLine).matches(rt.getStringValue(FILE_COL_NAME, rtLine-1))
					|| !rt.getStringValue(ROI_COL_NAME, rtLine).matches(rt.getStringValue(ROI_COL_NAME, rtLine-1))) {
				setPaths(directory+rt.getStringValue(FILE_COL_NAME, rtLine));
				openImageROIs(true, false);
			}
			if (verified && rt.getValue(VERIF_COL_NAME,rtLine) < 1) continue;
			float[] xPoints = getRM().getRoi((int)rt.getValue(ROI_COL_NAME, rtLine)-1).getFloatPolygon().xpoints;
			float[] yPoints = getRM().getRoi((int)rt.getValue(ROI_COL_NAME, rtLine)-1).getFloatPolygon().ypoints;
			double roiLength = rt.getValue(LENGTH_COL_NAME, rtLine);
			String pathXtCropFile = xtPathCrop+"XTcrop"+rt.getStringValue(ROI_COL_NAME, rtLine)+"_"+rt.getStringValue(FILE_COL_NAME, rtLine).replace(".tif", ".txt");
			float yStart = yPoints[0];
			float yStop = yPoints[yPoints.length-1];
			File xtCropFile = new File(pathXtCropFile);
			if (xtCropFile.exists()) {
				String cropParmsTxt = IJ.openAsString(pathXtCropFile);
				String[] cropParms = cropParmsTxt.split("-|\n");
				int yStartIdx = Integer.parseInt(cropParms[0]);
				int yStopIdx = Integer.parseInt(cropParms[1]);
				double distanceTravelled = 0;
				for (int iX = 0; iX < xPoints.length-1; iX++) {
					double distanceToAdd = Math.sqrt(Math.pow(xPoints[iX+1]-xPoints[iX], 2)+Math.pow(yPoints[iX+1]-yPoints[iX], 2));
					if (yStartIdx >= distanceTravelled && yStartIdx < distanceTravelled + distanceToAdd) {
						double frac = (yStartIdx - distanceTravelled) / distanceToAdd;
						yStart = (float)(yPoints[iX]*(1-frac) + yPoints[iX+1]*frac);
					}
					if (yStopIdx >= distanceTravelled && yStopIdx < distanceTravelled + distanceToAdd) {
						double frac = (yStopIdx - distanceTravelled) / distanceToAdd;
						yStop = (float)(yPoints[iX]*(1-frac) + yPoints[iX+1]*frac);
					}
					distanceTravelled += distanceToAdd;
				}
			}
			double dt = (Math.abs(camYOffset+yStop-rollingShutterStartLine)
					- Math.abs(camYOffset+yStart-rollingShutterStartLine)) * camLineTimeSec;
			double angleTan = Math.tan(Math.PI/180.0*rt.getValue(ANGLE_OK_COL_NAME, rtLine));
			
			double t = roiLength * angleTan * rt.getValue(XPX_MS_COL_NAME, rtLine) / 1000;
			double rsCoeff = t/(t+dt);
			double newVelocity = rt.getValue(YPX_UM_COL_NAME,rtLine)*1000/(rt.getValue(XPX_MS_COL_NAME,rtLine)*angleTan) * rsCoeff;
			//double newVelocity = (rt.getValue(VELOCITY_COL_NAME,rtLine) * Math.tan(Math.PI/180.0*rt.getValue(ANGLE_COL_NAME, rtLine)) / angleTan) / rsCoeff;
			rt.setValue(VELOCITY_OK_COL_NAME, rtLine, newVelocity);
			rt.setValue(RS_COEFF_COL_NAME, rtLine, rsCoeff);
		}
		rt.save(rtPathChoice);
		rt.show(velocityFileName);
	}
	
	private void drawVelocity() {
		if (emptyListAbort()) return;
		if (rtPath == null) return;
		
		String[] velocityFileNames = getExistingTableNames("All");
		if (velocityFileNames.length == 0) {
			new MessageDialog(new Frame(), "No result table", "No result table found. Aborting.");
			return;
		}
		
		GenericDialog gd = new GenericDialog("Draw velocity");
		gd.addChoice("Choose table:", velocityFileNames, velocityFileNames[0]);
		gd.addNumericField("Lower boundary (µm/s)", Prefs.get("RBCv.colorLowBoundary", 0));
		gd.addNumericField("Upper boundary (µm/s)", Prefs.get("RBCv.colorUpBoundary", 2000));
		gd.addChoice("Groups to display", DRAW_VELO_GROUP_CHOICES, Prefs.get("RBCv.colorUpBoundary",DRAW_VELO_GROUP_CHOICES[0]));
		gd.showDialog();
		if (gd.wasCanceled()) return;
		String velocityFileName = gd.getNextChoice();
		double boundaryLow = gd.getNextNumber();
		Prefs.set("RBCv.colorLowBoundary", boundaryLow);
		double boundaryUp = gd.getNextNumber();
		Prefs.set("RBCv.colorUpBoundary", boundaryUp);
		int drawGroup = (int)gd.getNextChoiceIndex();
		Prefs.set("RBCv.colorGroupChoice", DRAW_VELO_GROUP_CHOICES[drawGroup]);
		String tableName;
		
		String rtPathChoice;
		if (velocityFileName.matches(VELOCITY_FILENAME)) {
			rtPathChoice = rtPath;
			tableName = "_raw";
		} else if (velocityFileName.matches(VELOCITY_MANUAL_FILENAME)) {
			rtPathChoice = rtManualPath;
			tableName = "_manual";
		} else if (velocityFileName.matches(VELOCITY_AUTO_FILENAME)) {
			rtPathChoice = rtAutoPath;
			tableName = "_auto";
		} else if (velocityFileName.matches(VELOCITY_SLIDING_FILENAME)) {
			rtPathChoice = rtSlidingPath;
			tableName = "_sliding";
		} else { // if (velocityFileName.matches(VELOCITY_SLIDING_AUTO_FILENAME)) {
			rtPathChoice = rtSlidingAutoPath;
			tableName = "_sliding_auto";
		}
		try {
			rt = ResultsTable.open(rtPathChoice);
			rt.setPrecision(TABLE_PRECISION);
			rt.showRowNumbers(false);
		} catch (IOException e) {
			new MessageDialog(new Frame(), "IO Exception", "Could not open file\n"+rtPathChoice);
			return;
		}
		
		//rt.sort(PREPROCESS_COL_NAME);
		rt.sort(FILE_COL_NAME);
		rt.show(velocityFileName);
		boolean verified = (rt.getColumnIndex(VERIF_COL_NAME) != ResultsTable.COLUMN_NOT_FOUND);
		boolean hasOkCol = (rt.getColumnIndex(VELOCITY_OK_COL_NAME) != ResultsTable.COLUMN_NOT_FOUND);
		long startTime = System.currentTimeMillis();
		int[] rtEdgesTmp = new int[rt.size()+1];
		int rtEdgesIdx = 0;
		rtEdgesTmp[rtEdgesIdx++] = 0;
		int rtLine = 1;
		for (; rtLine < rt.size(); rtLine++) {
			if (!(rt.getStringValue(FILE_COL_NAME, rtLine-1).matches(rt.getStringValue(FILE_COL_NAME, rtLine))) ) {
					//|| !(rt.getStringValue(PREPROCESS_COL_NAME, rtLine-1).matches(rt.getStringValue(PREPROCESS_COL_NAME, rtLine)))) {
				rtEdgesTmp[rtEdgesIdx++] = rtLine;
			}
		}
		int[][] colorMap = {{47,2,119},{46,2,120},{45,2,122},{44,2,123},{42,2,124},{41,2,125},{40,2,126},{39,2,127},{37,2,128},{36,1,129},{35,1,131},{33,1,132},
				{32,1,133},{30,1,134},{29,1,135},{27,1,136},{26,1,137},{24,1,138},{23,1,139},{21,1,141},{19,1,142},{18,1,143},{16,1,144},{14,1,145},{12,1,146},
				{10,1,147},{9,0,148},{7,0,149},{5,0,151},{3,0,152},{1,0,153},{0,1,154},{0,3,155},{0,5,156},{0,7,157},{0,9,158},{0,11,160},{0,13,161},{0,15,162},
				{0,18,163},{1,20,164},{2,23,165},{3,26,166},{3,28,167},{4,31,168},{5,33,169},{6,36,169},{7,39,170},{8,41,171},{8,44,172},{9,47,173},{10,49,174},
				{11,52,175},{12,55,176},{13,58,177},{14,60,178},{15,63,179},{16,66,180},{16,69,181},{17,71,182},{18,74,183},{19,77,184},{20,80,185},{21,82,186},
				{22,85,187},{23,88,188},{24,91,189},{25,94,190},{26,96,191},{26,99,191},{26,102,191},{26,105,191},{26,108,192},{26,111,192},{26,114,192},{26,116,192},
				{26,119,192},{26,122,193},{26,125,193},{26,128,193},{26,131,193},{26,134,193},{26,137,194},{26,140,194},{26,143,194},{26,146,194},{26,149,194},
				{26,152,195},{26,155,195},{25,158,195},{25,161,195},{25,164,195},{25,167,196},{25,170,196},{25,173,196},{25,176,196},{25,179,196},{25,182,197},
				{25,185,197},{25,188,197},{25,191,197},{25,194,197},{25,198,198},{25,198,195},{25,198,192},{25,198,190},{25,198,187},{25,199,184},{25,199,181},
				{25,199,178},{25,199,176},{25,199,173},{25,200,170},{24,200,167},{24,200,165},{24,200,162},{24,200,159},{24,201,156},{24,201,153},{24,201,150},
				{24,201,148},{24,202,143},{24,202,138},{24,202,133},{24,203,129},{24,203,124},{24,203,119},{24,204,114},{24,204,109},{24,204,104},{23,205,99},
				{23,205,94},{23,205,89},{23,206,84},{23,206,79},{23,206,74},{23,207,69},{23,207,64},{23,207,59},{23,208,54},{23,208,49},{23,208,44},{23,209,39},
				{22,209,33},{22,209,28},{22,210,23},{27,210,22},{32,210,22},{37,211,22},{42,211,22},{48,211,22},{53,212,22},{58,212,22},{63,212,22},{69,213,22},
				{74,213,21},{79,213,21},{85,214,21},{90,214,21},{96,214,21},{101,215,21},{107,215,21},{112,215,21},{118,216,21},{123,216,21},{129,216,21},{134,217,21},
				{140,217,20},{145,217,20},{148,217,20},{151,217,20},{153,218,20},{156,218,20},{159,218,20},{161,218,20},{164,218,20},{166,218,20},{169,219,20},
				{172,219,20},{174,219,20},{177,219,20},{180,219,20},{182,219,20},{185,219,20},{188,220,20},{190,220,20},{193,220,20},{196,220,20},{198,220,19},
				{201,220,19},{204,221,19},{206,221,19},{210,221,19},{213,221,19},{216,221,19},{220,221,19},{222,220,19},{222,217,19},{222,214,19},{222,212,19},
				{222,209,19},{223,206,19},{223,203,19},{223,200,19},{223,197,19},{223,194,19},{223,191,19},{224,188,19},{224,185,18},{224,182,18},{224,178,18},
				{224,175,18},{225,172,18},{225,169,18},{225,166,18},{225,163,18},{225,160,18},{225,157,18},{226,154,18},{226,151,18},{226,148,18},{226,145,18},
				{226,141,18},{227,138,18},{227,135,18},{227,132,18},{227,129,17},{227,126,17},{228,122,17},{228,119,17},{228,116,17},{228,113,17},{228,110,17},
				{228,106,17},{229,103,17},{229,100,17},{229,97,17},{229,94,17},{229,90,17},{230,87,17},{230,84,17},{230,80,17},{230,77,17},{230,74,17},{230,71,16},
				{231,67,16},{231,64,16},{231,61,16},{231,57,16},{231,54,16},{232,51,16},{232,47,16},{232,44,16},{232,41,16},{232,37,16},{233,34,16}};
		rtEdgesTmp[rtEdgesIdx++] = rt.size();
		int[] rtEdges = new int[rtEdgesIdx];
		System.arraycopy(rtEdgesTmp, 0, rtEdges, 0, rtEdgesIdx);
		ImagePlus velocityImp = new ImagePlus();
		for (int iEdge = 0; iEdge < rtEdges.length-1; iEdge++) {
			File f = new File(directory+rt.getStringValue(FILE_COL_NAME, rtEdges[iEdge]));
			String name = f.getName();
			name = name.substring(0, name.lastIndexOf("."))+tableName+".png";
			String velocityDirectory = f.getParentFile().getPath()+File.separator+"VelocityArrows"+File.separator;
			f = new File(velocityDirectory);
			if (!f.exists()) {
				f.mkdir();
			}
			if (iEdge == 0 || !(rt.getStringValue(FILE_COL_NAME, rtEdges[iEdge-1]).matches(rt.getStringValue(FILE_COL_NAME, rtEdges[iEdge])))) {
				if (oriImp != null) oriImp.close();
				if (velocityImp != null) velocityImp.close();
				oriImp = IJ.openImage(directory+rt.getStringValue(FILE_COL_NAME, rtEdges[iEdge]));
				setPaths(directory+rt.getStringValue(FILE_COL_NAME, rtEdges[iEdge]));
				velocityImp = ZProjector.run(oriImp,"avg");
				if (velocityImp.getType() != ImagePlus.GRAY32) {
					velocityImp.setProcessor(velocityImp.getProcessor().convertToFloat());
				}
				oriImp.close();
				Rectangle boundingBox = getNotNullBoundingBox(velocityImp);
				velocityImp.setRoi(boundingBox);
				ImagePlus velocityHiImp = velocityImp.crop();
				ImagePlus velocityLoImp = velocityImp.crop();
				
				double gaussHi = 10;
				double gaussLo = 2;
				IJ.run(velocityHiImp, "Gaussian Blur...", "sigma="+gaussHi);
				IJ.run(velocityLoImp, "Gaussian Blur...", "sigma="+gaussLo);
				ImageCalculator.run(velocityLoImp, velocityHiImp, "subtract");
				velocityLoImp.copy();
				velocityImp.paste();
				velocityImp.deleteRoi();
				velocityImp.getProcessor().resetMinAndMax();
				IJ.run(velocityImp, "Enhance Contrast", "saturated=0.35");
				velocityHiImp.close();
				velocityLoImp.close();
			}

			IJ.run(velocityImp, "Remove Overlay", "");
			openImageROIs(true, false);

			int nROI = getRM().getCount();
			int nROITable = rtEdges[iEdge+1] - rtEdges[iEdge];
			int nROITimes = nROITable / nROI;
			if (nROITable % nROI != 0) {
				IJ.log("Different number of ROIs in ROI Manager ("+nROI+") and in Table ("+nROITable+") for File "+
						rt.getStringValue(FILE_COL_NAME, rtEdges[iEdge])+", skipping it.");
				continue;
			}
			IJ.log("Image "+(iEdge+1)+"/"+(rtEdges.length-1)+", "+nROI+" ROIs. "+((System.currentTimeMillis()-startTime)/1000.0)+"s elapsed.");
			for (int iROI = 0; iROI < nROI; iROI++) {
				double verifValue = 0;
				if (verified) verifValue = rt.getValue(VERIF_COL_NAME, rtEdges[iEdge]+iROI);
				double velocityValue = 0;
				if (nROITimes > 1) {
					double nEffectiveTimes = 0;
					for (int iTimes = 0; iTimes < nROITimes; iTimes++) {
						int iTimesIndex = rtEdges[iEdge]+iROI*nROITimes+iTimes;
						if ((verified && rt.getValue(VERIF_COL_NAME, iTimesIndex) > 0) || !verified) {
							velocityValue += hasOkCol ? rt.getValue(VELOCITY_OK_COL_NAME, iTimesIndex) : rt.getValue(VELOCITY_COL_NAME, iTimesIndex);
							nEffectiveTimes += 1;
							if (verified) {
								verifValue = Math.max(verifValue, rt.getValue(VERIF_COL_NAME, iTimesIndex));
							}
						}
					}
					if (nEffectiveTimes > 0) {
						velocityValue /= nEffectiveTimes;
					}
				} else {
					if (hasOkCol) {
						velocityValue = rt.getValue(VELOCITY_OK_COL_NAME, rtEdges[iEdge]+iROI);
					} else {
						velocityValue = rt.getValue(VELOCITY_COL_NAME, rtEdges[iEdge]+iROI);
					}
				}
				if (verifValue > -1) {
					FloatPolygon fp = toPolyline(getRM().getRoi(iROI).getFloatPolygon(), 20);
					fp = toArrowPolyline(fp, velocityValue < 0);
					PolygonRoi pr = new PolygonRoi(fp, Roi.POLYLINE);
					int groupVal = (int)Math.pow(2, rt.getValue(GROUP_COL_NAME, rtEdges[iEdge]+iROI));
					if ((groupVal & DRAW_VELO_GROUP_CHOICES_VAL[drawGroup]) > 0) {
						//pr.setStrokeColor(Color.YELLOW);
						int indexColor = (int)((colorMap.length - 1) * (Math.abs(velocityValue) - boundaryLow) / (boundaryUp - boundaryLow));
						if (indexColor < 0) indexColor = 0;
						if (indexColor > (colorMap.length - 1)) indexColor = (colorMap.length - 1);
						pr.setStrokeColor(new Color(colorMap[indexColor][0],colorMap[indexColor][1],colorMap[indexColor][2]));
						//pr.setName(""+(iROI+1)+"-"+Math.round(Math.abs(velocityValue)));
						pr.setStrokeWidth(getRM().getRoi(iROI).getStrokeWidth());
						velocityImp.setRoi(pr);
						IJ.run(velocityImp, "Add Selection...", "");
					}
				}
			}
			
			IJ.run(velocityImp, "Labels...", "color=green font=16 show use bold");
			ImagePlus flattenImp = velocityImp.flatten();
			IJ.saveAs(flattenImp, "png", velocityDirectory+name);
			flattenImp.close();
		}
		velocityImp.close();
		IJ.log("Draw velocity done. "+((System.currentTimeMillis()-startTime)/1000.0)+"s elapsed.");
	}
	
	public static FloatPolygon toPolyline(FloatPolygon fp, float spacing) {
		// change current ROI to a polyline with equal segments of size approximately "spacing"
		float[] xIn = fp.xpoints;
		float[] yIn = fp.ypoints;
		int nPointIn = xIn.length;
		float[] cumulDistance = new float[nPointIn];
		for (int iPoint = 1; iPoint < nPointIn; iPoint++) {
			cumulDistance[iPoint] = cumulDistance[iPoint-1] + (float)Math.sqrt(Math.pow(xIn[iPoint]-xIn[iPoint-1],2)+Math.pow(yIn[iPoint]-yIn[iPoint-1],2));
		}
		
		int nPointOut = (int)Math.round(cumulDistance[nPointIn-1]/spacing);
		spacing = cumulDistance[nPointIn-1]/nPointOut;
		float[] xOut = new float[nPointOut+1];
		float[] yOut = new float[nPointOut+1];
		xOut[0] = xIn[0];
		yOut[0] = yIn[0];
		
		int idx = 1;
		for (int iPoint = 0; iPoint < nPointOut; iPoint++) {
			float pos = iPoint*spacing;
			while (cumulDistance[idx] < pos) {
				idx++;
			}
			float coeff1 = (cumulDistance[idx]-pos) / (cumulDistance[idx]-cumulDistance[idx-1]);
			xOut[iPoint] = coeff1 * xIn[idx-1] + (1-coeff1) * xIn[idx];
			yOut[iPoint] = coeff1 * yIn[idx-1] + (1-coeff1) * yIn[idx];
		}
		xOut[nPointOut] = xIn[nPointIn-1];
		yOut[nPointOut] = yIn[nPointIn-1];
		return new FloatPolygon(xOut,yOut);
	}
	
	public static FloatPolygon toArrowPolyline(FloatPolygon fp, boolean reverseOrder) {
		float[] xIn = fp.xpoints;
		float[] yIn = fp.ypoints;
		float[] xOut = new float[fp.npoints+3];
		float[] yOut = new float[fp.npoints+3];
		if (reverseOrder) {
			for (int i = 0; i < fp.npoints; i++) {
				xOut[i] = xIn[fp.npoints - 1 - i];
				yOut[i] = yIn[fp.npoints - 1 - i];
			}
		} else {
			for (int i = 0; i < fp.npoints; i++) {
				xOut[i] = xIn[i];
				yOut[i] = yIn[i];
			}
		}
		xOut[fp.npoints+2] = xOut[fp.npoints-1];
		yOut[fp.npoints+2] = yOut[fp.npoints-1];
		
		float dx = xOut[fp.npoints-1]-xOut[fp.npoints-2];
		float dy = yOut[fp.npoints-1]-yOut[fp.npoints-2];
		float dxy = (float)Math.sqrt(dx*dx+dy*dy);
		dx = dx/dxy;
		dy = dy/dxy;
		float a = 12;
		xOut[fp.npoints] = a*(-dx+dy)+xOut[fp.npoints-1];
		xOut[fp.npoints+1] = a*(-dx-dy)+xOut[fp.npoints-1];
		yOut[fp.npoints] = a*(-dy-dx)+yOut[fp.npoints-1];
		yOut[fp.npoints+1] = a*(-dy+dx)+yOut[fp.npoints-1];
		return new FloatPolygon(xOut, yOut);
	}
	
	private Rectangle getNotNullBoundingBox(ImagePlus imp) {
		int w = imp.getWidth();
		int h = imp.getHeight();
		Rectangle res = new Rectangle();
		for (int iW = 0; iW < w; iW++) {
			imp.setRoi(new Rectangle(iW,0,1,h));
			ImageStatistics is = imp.getStatistics(ImageStatistics.MIN_MAX);
			if (is.max > 0) {
				res.x = iW;
				break;
			}
		}
		for (int iW = w-1; iW >= 0; iW--) {
			imp.setRoi(new Rectangle(iW,0,1,h));
			ImageStatistics is = imp.getStatistics(ImageStatistics.MIN_MAX);
			if (is.max > 0) {
				res.width = Math.max(iW-res.x+1, 1);
				break;
			}
		}
		for (int iH = 0; iH < h; iH++) {
			imp.setRoi(new Rectangle(0,iH,w,1));
			ImageStatistics is = imp.getStatistics(ImageStatistics.MIN_MAX);
			if (is.max > 0) {
				res.y = iH;
				break;
			}
		}
		for (int iH = h-1; iH >= 0; iH--) {
			imp.setRoi(new Rectangle(0,iH,w,1));
			ImageStatistics is = imp.getStatistics(ImageStatistics.MIN_MAX);
			if (is.max > 0) {
				res.height = Math.max(iH-res.y+1, 1);
				break;
			}
		}
		imp.resetRoi();
		return res;
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Object b = e.getSource();
		if (b == selectDirectoryBtn) {
			updateDirectory();
		} else if (b == stackSuffixValue && stackSuffixChanged) {
			updateStacksList();
			stackSuffixChanged = false;
		//} else if (b == shiftSuffixValue && shiftSuffixChanged) {
		//	updateStacksList();
		//	shiftSuffixChanged = false;
		} else if (b == openStackBtn) {
			openStack(false, true);
		} else if (b == saveROIsBtn) {
			saveROIs();
		} else if (b == exportROIsBtn) {
			exportROIsToOtherStacks();
		} else if (b == xtGenBtn) {
			Thread t1 = new Thread(new Runnable() {
				public void run() {
					xtGen();
				}
			});
			t1.start();
		} else if (b == measureVelocityBtn) {
			Thread t1 = new Thread(new Runnable() {
				public void run() {
					measureVelocity(false);
				}
			});
			t1.start();
		} else if (b == measureVelocitySlidingBtn) {
			Thread t1 = new Thread(new Runnable() {
				public void run() {
					measureVelocity(true);
				}
			});
			t1.start();
		} else if (b == manualVerifBtn) {
			manualVerification();
		} else if (b == evaluateVerifBtn) {
			evaluateVerification();
		} else if (b == automaticVerifBtn) {
			automaticVerification();
		} else if (b == rollingShutterBtn) {
			rollingShutterAdjust();
		} else if (b == velocityLimitBtn) {
			velocityLimit();
		} else if (b == drawVelocityBtn) {
			Thread t1 = new Thread(new Runnable() {
				public void run() {
					drawVelocity();
				}
			});
			t1.start();
		} else if (b == updateAngleBtn) {
			updateAngle();
		} else if (b == discardAngleBtn) {
			discardAngle();
		} else if (b == openStackVerifBtn) {
			openStackVerif();
		} else if (b == openXTRawVerifBtn) {
			openXTRawVerif();
		} else if (b == validateROIBtn) {
			validateROI();
		} else if (b == unvalidateROIBtn) {
			unvalidateROI();
		}
	}
	@Override
	public void mouseClicked(MouseEvent e) {}
	@Override
	public void mousePressed(MouseEvent e) {}
	@Override
	public void mouseReleased(MouseEvent e) {}
	@Override
	public void mouseEntered(MouseEvent e) {}
	@Override
	public void mouseExited(MouseEvent e) {}
	@Override
	public void windowOpened(WindowEvent e) {}
	@Override
	public void windowClosing(WindowEvent e) {
		Object b = e.getSource();
		if (b == frame) {
			closeImps(true);
			closeVerif();
			closeFrame();
		} else if (b == verif) {
			closeImps(false);
			closeVerif();
		}
	}
	@Override
	public void windowClosed(WindowEvent e) {}
	@Override
	public void windowIconified(WindowEvent e) {}
	@Override
	public void windowDeiconified(WindowEvent e) {}
	@Override
	public void windowActivated(WindowEvent e) {}
	@Override
	public void windowDeactivated(WindowEvent e) {}
	@Override
	public void textValueChanged(TextEvent e) {
		Object b = e.getSource();
		if (b == stackSuffixValue) {
			stackSuffixChanged = true;
		//} else if (b == shiftSuffixValue) {
		//	shiftSuffixChanged = true;
		}
	}
	@Override
	public void itemStateChanged(ItemEvent e) {
		if (e.getSource() == depthValue) {
			updateStacksList();
		}
	}
	@Override
	public void focusGained(FocusEvent e) {}
	@Override
	public void focusLost(FocusEvent e) {
		Object b = e.getSource();
		if (b == stackSuffixValue && stackSuffixChanged) {
			updateStacksList();
			stackSuffixChanged = false;
		//} else if (b == shiftSuffixValue && shiftSuffixChanged) {
		//	updateStacksList();
		//	shiftSuffixChanged = false;
		}
	}
	
	private void closeFrame() {
		if (frame != null) {
			framePosX = frame.getLocation().getX();
			framePosY = frame.getLocation().getY();
			ij.Prefs.set("RBCv.framePosX", framePosX);
			ij.Prefs.set("RBCv.framePosY", framePosY);
			frame.dispose();
			frame = null;
			Hashtable shortcuts = Menus.getShortcuts();
			for (int iKeys = 0; iKeys < setROIsShortcutsKeys.length; iKeys++) {
				int keyValue = Menus.convertShortcutToCode(setROIsShortcutsKeys[iKeys]);
				if (setROIsShortcutsCode[iKeys] != null) {
					shortcuts.put(keyValue,setROIsShortcutsCode[iKeys]);
				} else {
					shortcuts.remove(keyValue);
				}
			}
			for (int iKeys = 0; iKeys < groupColorShortcutsKeys.length; iKeys++) {
				int keyValue = Menus.convertShortcutToCode(groupColorShortcutsKeys[iKeys]);
				if (groupColorShortcutsCode[iKeys] != null) {
					shortcuts.put(keyValue,groupColorShortcutsCode[iKeys]);
				} else {
					shortcuts.remove(keyValue);
				}
			}
		}
	}
	
	private void closeVerif() {
		if (verif != null) {
			verifPosX = verif.getLocation().getX();
			verifPosY = verif.getLocation().getY();
			ij.Prefs.set("RBCv.verifPosX", verifPosX);
			ij.Prefs.set("RBCv.verifPosY", verifPosY);
			verif.dispose();
			verif = null;
			Hashtable shortcuts = Menus.getShortcuts();
			for (int iKeys = 0; iKeys < verifShortcutsKeys.length; iKeys++) {
				int keyValue = Menus.convertShortcutToCode(verifShortcutsKeys[iKeys]);
				if (verifShortcutsCode[iKeys] != null) {
					shortcuts.put(keyValue,verifShortcutsCode[iKeys]);
				} else {
					shortcuts.remove(keyValue);
				}
			}
		}
	}
	
	private void closeImps(boolean everything) {
		if (xtImp != null) {xtImp.close(); xtImp = null;}
		if (verifImp != null) {verifImp.close(); verifImp = null;}
		if (everything) {
			if (oriImp != null) {oriImp.close(); oriImp = null;}
			if (allStacksImp != null) {allStacksImp.close(); allStacksImp = null;}
			if (xtPreImp != null) {
				for (int iXT = 0; iXT < xtPreImp.length; iXT++) {
					if (xtPreImp[iXT] != null) xtPreImp[iXT].close();
				}
				xtPreImp = null;
			}
		}
		getRM().close();
	}
	
	private boolean openImageROIs(boolean reset, boolean showAll) {
		if (reset) {
			getRM().reset();
		}
		File f = new File(oriROIsPath);
		if (f.exists()) {
			getRM().runCommand("Open", oriROIsPath);
			if (showAll) {
				getRM().runCommand(oriImp, "Show All with Labels");
			}
			return true;
		} else {
			return false;
		}
	}
	
	private void radonTransformDialog() {
		ImagePlus imp = IJ.getImage();
		
		GenericDialog gd = new GenericDialog("Radon transform parameters");
		gd.addNumericField("Min angle (included)", 0, 1, 8, "");
		gd.addNumericField("Step angle", 1, 1, 8, "");
		gd.addNumericField("Max angle (excluded)", 180, 1, 8, "");
		gd.addNumericField("Subdivisions", 100, 1, 8, "");
		gd.showDialog();
		double minAngle = gd.getNextNumber();
		double stepAngle = gd.getNextNumber();
		double maxAngle = gd.getNextNumber();
		int subdivisions = (int)gd.getNextNumber();
		if (stepAngle == 0) stepAngle = 1;
		int nAngles = (int)Math.ceil((maxAngle-minAngle) / stepAngle);
		double[] angles = new double[nAngles];
		for (int i = 0; i < nAngles; i++) {
			angles[i] = minAngle + i * stepAngle;
		}
		ImagePlus imp2 = radonTransform(imp, angles, subdivisions);
		IJ.resetMinAndMax(imp2);
		imp2.show();
	}
	
	public static ImagePlus radonTransform(ImagePlus imp) {
		double[] angles = new double[180];
		for (int i = 0; i < 180; i++) {
			angles[i] = i;
		}
		return radonTransform(imp, angles, 100);
	}
	
	public static ImagePlus radonTransform(ImagePlus imp, double[] angles) {
		return radonTransform(imp, angles, 100);
	}

	public static ImagePlus radonTransform(ImagePlus imp, double[] angles, int subdivisions) {
		ImageProcessor ip = imp.getChannelProcessor();
		int radonWidth = angles.length;
		int imageWidth = imp.getWidth();
		int imageHeight = imp.getHeight();
		int radonHeight = (int)Math.ceil(Math.sqrt(imageHeight*imageHeight + imageWidth*imageWidth)) + 4;
		double radonHeightHalf = radonHeight / 2;
		double[] cosvalCoeff = new double[imageWidth]; // x-axis distance from pixel to image center
		double[] sinvalCoeff = new double[imageHeight]; // y-axis distance from pixel to image center
		ImagePlus imp2 = IJ.createImage("radon_transformed", "32-bit black", radonWidth, radonHeight, 1);
		ImageProcessor ip2 = imp2.getProcessor();
		int index = 0;
		for (double value = -(imageWidth-1)/2; value <= (imageWidth-1)/2; value++) {
			cosvalCoeff[index++] = value;
		}
		index = 0;
		for (double value = -(imageHeight-1)/2; value <= (imageHeight-1)/2; value++) {
			sinvalCoeff[index++] = value;
		}
		
		IntStream.range(0, radonWidth).parallel().forEach(iAngle -> {
			double[] radonColumnSum = new double[radonHeight];
			// 180-angle to make it angle from vertical downward axis, counterclockwise
			double angle = (180-angles[iAngle])*Math.PI/180;
			double cosval = Math.cos(angle);
			double sinval = Math.sin(angle);
			
			double[] pxVerticesProj = new double[4];
			pxVerticesProj[0] = (-0.5 * cosval -0.5 * sinval + 2.0) * subdivisions;
			pxVerticesProj[1] = (-0.5 * cosval +0.5 * sinval + 2.0) * subdivisions;
			pxVerticesProj[2] = (+0.5 * cosval -0.5 * sinval + 2.0) * subdivisions;
			pxVerticesProj[3] = (+0.5 * cosval +0.5 * sinval + 2.0) * subdivisions;
			Arrays.sort(pxVerticesProj);
			double[] bins = new double[4 * subdivisions + 1];
			if (pxVerticesProj[1] != pxVerticesProj[0]) {
				for (int idx = (int)Math.ceil(pxVerticesProj[0]); idx <= (int)pxVerticesProj[1]; idx++) {
					bins[idx] = (idx-pxVerticesProj[0])/(pxVerticesProj[1]-pxVerticesProj[0]);
				}
			}
			for (int idx = (int)Math.ceil(pxVerticesProj[1]); idx <= (int)pxVerticesProj[2]; idx++) {
				bins[idx] = 1;
			}
			if (pxVerticesProj[3] != pxVerticesProj[2]) {
				for (int idx = (int)Math.ceil(pxVerticesProj[2]); idx <= (int)pxVerticesProj[3]; idx++) {
					bins[idx] = (pxVerticesProj[3]-idx)/(pxVerticesProj[3]-pxVerticesProj[2]);
				}
			}
			// normalize to a sum of 1
			double binsSum = 0;
			for (int idx = (int)pxVerticesProj[0]; idx <= (int)pxVerticesProj[3]; idx++) {
				binsSum += bins[idx];
			}
			for (int idx = (int)pxVerticesProj[0]; idx <= (int)pxVerticesProj[3]; idx++) {
				bins[idx] /= binsSum;
			}
			// cumulative sum
			for (int idx = 1; idx < bins.length; idx++) {
				bins[idx] += bins[idx-1];
			}
			// derivation with step "subdivisions"
			double[] deltaBins = new double[3*subdivisions+1];
			for (int i = 0; i < deltaBins.length; i++) {
				deltaBins[i] = bins[i+subdivisions] - bins[i];
			}
			for (int x = 0; x < imageWidth; x++) {
				double cosvalXCoeff = cosval * cosvalCoeff[x] + radonHeightHalf;
				for (int y = 0; y < imageHeight; y++) {
					double pos = cosvalXCoeff + sinval * sinvalCoeff[y];
					int pos1 = (int)pos;
					int c1 = (int)(subdivisions*(pos-pos1)); // eventually add rounding here
					double pxVal = ip.getPixelValue(x, y);
					radonColumnSum[pos1+1] += pxVal * deltaBins[c1];
					radonColumnSum[pos1] += pxVal * deltaBins[c1+subdivisions];
					radonColumnSum[pos1-1] += pxVal * deltaBins[c1+2*subdivisions];
				}
			}
			for (int iCol = 0; iCol < radonColumnSum.length; iCol++) {
				ip2.putPixelValue(iAngle,iCol, radonColumnSum[iCol]);
			}
	    });
		return imp2;
	}
	
	private RoiManager getRM() {
		return RoiManager.getRoiManager();
	}
	
	/*
	private static String getCommand(String cmdVal) {
		String res = IJ.runMacro("str = getArgument(); List.setCommands; List.toArrays(keys, values);"+
		"for (i = 0; i < values.length; i++) {if (endsWith(values[i], str)) {return keys[i];}}", cmdVal);
		if (res == null) {
			res = cmdVal+" plugin not found.\nInstall the plugin and restart the application to use it.";
		}
		return res;
	}
	*/
}



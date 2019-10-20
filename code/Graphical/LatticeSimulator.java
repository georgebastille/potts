
/*
 *  An implementation of the q-state Potts model
 *  Uses a square lattice with periodic boundary conditions.
 *  Metropolis & Wolff Dynamic Algorithm implemented.
 *  NB temp is in units of Kelvin^2/Joule and J = kB = 1
 *
 *  Richard Hanes		17/12/04
 */
import java.io.*;
import java.awt.*;
import java.awt.event.*;

/**
 *  Contains all the information and methods to manipulate a Lattice
 *
 * @author    Richard Hanes
 */
class Lattice {
	boolean running = true;
	boolean metropolis = false;
	boolean sizeChanged = false, statesChanged = false, resetMe = false;
	int size, states, systemEnergy, newSize, newStates;
	int[] magnetisation;
	/**
	 *  The 2d array which holds the current state of each point in the lattice
	 */
	public int latticeArray[][];
	boolean[][] cluster;

	double exponentials[];
	double temp, tempGradient, addProb;


	/**
	 *Constructor for the Lattice object
	 *
	 * @param  size          Linear Size of the Lattice
	 * @param  states        No of Spin states in the simulation
	 * @param  temp          Initial temperature of the Lattice
	 * @param  tempGradient  Initial Heating/Cooling Gradient
	 */
	public Lattice(int size, int states, double temp, double tempGradient) {
		this.size = size;
		this.states = states;
		this.temp = temp;
		this.tempGradient = tempGradient;
		this.latticeArray = new int[size][size];
		//Used to store Metropolis transition Probabilities
		this.exponentials = new double[4];
		// Stores number of points in each state
		this.magnetisation = new int[states];
		// Wolff cluster add probability
		this.addProb = (1.0 - Math.exp(-1.0 / temp));
		this.cluster = new boolean[size][size];
		resetSimulation();
	}


	/**
	 *  Re-Read Lattice Parameters and Reset the lattice
	 */
	public void resetSimulation() {
		resetMagnetArray();
		if (metropolis == false) {
			resetClusterArray();
		}
		calculateExponentials();
		addProb = (1.0 - Math.exp(-1.0 / temp));
		randomiseLattice();
		calculateSystemMagnetisation();
		systemEnergy = calculateTotalSystemEnergy();
	}


	/**
	 *  Reset the Boolean Cluster Array (Wolff Algorithm)
	 */
	public void resetClusterArray() {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				cluster[i][j] = false;
			}
		}
	}


	/**
	 *  Reset the Magnetisation Array, where the number of points in each state is stored
	 */
	public void resetMagnetArray() {
		for (int i = 0; i < states; i++) {
			magnetisation[i] = 0;
		}
	}


	/**
	 *  Count the number of points in each state
	 */
	public void calculateSystemMagnetisation() {
		int currentPointState;

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				currentPointState = latticeArray[i][j];
				if (currentPointState <= states) {
					magnetisation[currentPointState - 1]++;
				}
			}
		}
	}


	/**
	 *  Calculate the flip probabilities (Metropolis Algorithm)
	 */
	public void calculateExponentials() {
		for (int i = 1; i <= 4; i++) {
			exponentials[i - 1] = Math.exp(-i / temp);
		}
	}


	/**
	 *  Set each point on the lattice to a random state (Hot Start)
	 */
	public void randomiseLattice() {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// Load a random state into each lattice point
				this.latticeArray[i][j] = (int) (Math.random() * states) + 1;
			}
		}
	}


	/**
	 *  Go through the entire lattice and calculate the energy from each point
	 *
	 * @return    Total System Energy
	 */
	public int calculateTotalSystemEnergy() {
		int energy = 0;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				energy += calculateEnergy(i, j);
			}
		}
		return energy;
	}


	/**
	 *  Calculate the Energy of a single lattice point
	 *
	 * @param  x  Coordinate of Lattice Point
	 * @param  y  Coordinate of Lattice Point
	 * @return    The energy of the point at (X,Y)
	 */
	public int calculateEnergy(int x, int y) {
		int[] nearestNeighboursArray = new int[4];
		int targetCellState = latticeArray[x][y];
		int energy = 0;

		//Periodic Boundary Conditions
		if (y == 0) {
			nearestNeighboursArray[0] = latticeArray[x][size - 1];
		} else {
			nearestNeighboursArray[0] = latticeArray[x][y - 1];
		}

		if (x == 0) {
			nearestNeighboursArray[1] = latticeArray[size - 1][y];
		} else {
			nearestNeighboursArray[1] = latticeArray[x - 1][y];
		}

		if (y == (size - 1)) {
			nearestNeighboursArray[2] = latticeArray[x][0];
		} else {
			nearestNeighboursArray[2] = latticeArray[x][y + 1];
		}

		if (x == (size - 1)) {
			nearestNeighboursArray[3] = latticeArray[0][y];
		} else {
			nearestNeighboursArray[3] = latticeArray[x + 1][y];
		}

		// This does the work of the Kronecker d-function
		for (int i = 0; i < 4; i++) {
			if (nearestNeighboursArray[i] == targetCellState) {
				energy++;
			}
		}

		return (-1 * energy);
	}



	/**
	 *  Run the Metropolis Algorithm
	 *
	 * @return    The Lattice point which was selected, for Drawing the Canvas
	 */
	public LatticePoint runMetropolis() {
		int oldEnergy;
		int newEnergy;
		int energyDifference;
		int systemEnergyDifference;
		LatticePoint chosen;

		// Pick a random lattice point to Flip
		int xToFlip = (int) (Math.random() * size);
		int yToFlip = (int) (Math.random() * size);

		// Remember the old state in case we are gonna flip back
		int oldState = latticeArray[xToFlip][yToFlip];

		// Pick a different random state to flip it to
		int flipToState = (int) (Math.random() * states) + 1;
		while (flipToState == oldState) {
			flipToState = (int) (Math.random() * states) + 1;
		}

		chosen = new LatticePoint(xToFlip, yToFlip, flipToState, true);

		// Calculate old state Energy
		oldEnergy = calculateEnergy(xToFlip, yToFlip);

		// Flip the state
		latticeArray[xToFlip][yToFlip] = flipToState;

		// Calculate the Energy of the new flipped state
		newEnergy = calculateEnergy(xToFlip, yToFlip);

		// This variable can take on values from -noOfNeighbours to noOfNeighbours
		energyDifference = newEnergy - oldEnergy;

		// The difference in energy of the system is twice the flipped cell difference
		systemEnergyDifference = 2 * energyDifference;

		// If the system is more stable (lower energy) continue
		if (energyDifference < 1) {
			systemEnergy += systemEnergyDifference;

			if (oldState <= states) {
				magnetisation[oldState - 1]--;
			}
			magnetisation[flipToState - 1]++;
			return chosen;
		} else {
			// Otherwise flip with a probability based on temperature
			double flipProb = Math.random();
			if (tempGradient != 0) {
				calculateExponentials();
			}

			if (exponentials[energyDifference - 1] > flipProb) {
				systemEnergy += systemEnergyDifference;
				if (oldState <= states) {
					magnetisation[oldState - 1]--;
				}
				magnetisation[flipToState - 1]++;
				return chosen;
			} else {
				latticeArray[xToFlip][yToFlip] = oldState;
				chosen.setState(oldState);
				chosen.setChanged(false);
				return chosen;
			}
		}
	}


	/**
	 *  Run the Wolff Algorithm
	 *
	 * @return    The Cluster which was flipped, for drawing the canvas
	 */
	public LatticeCluster runWolff() {
		addProb = (1.0 - Math.exp(-1.0 / temp));

		resetClusterArray();

		int xToFlip = (int) (Math.random() * size);
		int yToFlip = (int) (Math.random() * size);

		int oldState = latticeArray[xToFlip][yToFlip];

		int newState = (int) (Math.random() * states) + 1;
		while (newState == oldState) {
			newState = (int) (Math.random() * states) + 1;
		}

		growCluster(xToFlip, yToFlip, oldState, newState);
		LatticeCluster safe = new LatticeCluster(oldState, newState, cluster);
		return safe;
	}


	/**
	 *  Flip the selected point and try to add its neighbours to the cluster
	 *
	 * @param  x         Coordinate of the point to flip
	 * @param  y         Coordinate of the point to flip
	 * @param  oldState  The initial state of the point to flip
	 * @param  newState  The new state of the point to flip
	 */
	private void growCluster(int x, int y, int oldState, int newState) {
		cluster[x][y] = true;

		int oldEnergy = calculateEnergy(x, y);

		latticeArray[x][y] = newState;

		int newEnergy = calculateEnergy(x, y);

		int energyDiff = newEnergy - oldEnergy;
		systemEnergy += 2 * energyDiff;

		int xPrev = (x == 0) ? size - 1 : x - 1;
		int xNext = (x == size - 1) ? 0 : x + 1;
		int yPrev = (y == 0) ? size - 1 : y - 1;
		int yNext = (y == size - 1) ? 0 : y + 1;

		if (!cluster[xPrev][y]) {
			tryAdd(xPrev, y, oldState, newState);
		}
		if (!cluster[xNext][y]) {
			tryAdd(xNext, y, oldState, newState);
		}
		if (!cluster[x][yPrev]) {
			tryAdd(x, yPrev, oldState, newState);
		}
		if (!cluster[x][yNext]) {
			tryAdd(x, yNext, oldState, newState);
		}
	}


	/**
	 *  Try to add a Lattice point to the cluster
	 *
	 * @param  x         Coordinate of the point to flip
	 * @param  y         Coordinate of the point to flip
	 * @param  oldState  The initial state of the point to flip
	 * @param  newState  The new state of the point to flip
	 */
	private void tryAdd(int x, int y, int oldState, int newState) {
		if (latticeArray[x][y] == oldState) {
			if (Math.random() < addProb) {
				growCluster(x, y, oldState, newState);
			}
		}
	}


	/**
	 *  Change the Lattice Parameters
	 *
	 * @param  size    The new size of the Lattice
	 * @param  states  The new number of states any given point can take
	 */
	public void changeLattice(int size, int states) {
		this.size = size;
		this.states = states;
		this.latticeArray = new int[size][size];
		this.exponentials = new double[4];
		this.magnetisation = new int[states];
	}


	/**
	 *  Gets the size of the Lattice object
	 *
	 * @return    The size value
	 */
	public int getSize() {
		return this.size;
	}


	/**
	 *  Gets the number of states available to each spin in the Lattice
	 *
	 * @return    The no of states
	 */
	public int getStates() {
		return this.states;
	}


	/**
	 *  Change the size of the lattice while the simulation is running
	 *  If the size if decreasing, discard the hidden points
	 *  If the size is increasing, randomise the extra points
	 */
	public void changeLatticeSize() {
		int[][] newLattice = new int[newSize][newSize];
		cluster = new boolean[newSize][newSize];

		if (newSize < size) {
			for (int i = 0; i < newSize; i++) {
				for (int j = 0; j < newSize; j++) {
					newLattice[i][j] = latticeArray[i][j];
				}
			}
			latticeArray = newLattice;
			size = newSize;
			systemEnergy = calculateTotalSystemEnergy();
		} else {
			int[][] oldLattice = latticeArray;
			latticeArray = newLattice;
			int oldSize = size;
			size = newSize;
			randomiseLattice();
			for (int i = 0; i < oldSize; i++) {
				for (int j = 0; j < oldSize; j++) {
					latticeArray[i][j] = oldLattice[i][j];
				}
			}
			systemEnergy = calculateTotalSystemEnergy();
		}
	}


	/**
	 *  Change the number of states in the lattice and recalculate the no of points in each state
	 */
	public void changeNoOfStates() {
		int[] newMagArray = new int[newStates];
		magnetisation = newMagArray;
		states = newStates;
		calculateSystemMagnetisation();

	}
}

/**
 *  A small object to hold the coordinates, state and whether a Lattice point was flipped
 *
 * @author    Richard Hanes
 */
class LatticePoint {


	int x, y, state;
	boolean changed;


	/**
	 *Constructor for the LatticePoint object
	 *Used for for the Canvas to draw Metropolis flipped spins
	 *
	 * @param  x        Coordinate of Lattice Point in the Lattice
	 * @param  y        Coordinate of Lattice Point in the Lattice
	 * @param  state    State of the point
	 * @param  changed  Whether or not the point has changed
	 */
	public LatticePoint(int x, int y, int state, boolean changed) {
		this.x = x;
		this.y = y;
		this.state = state;
		this.changed = changed;
	}


	/**
	 *  Gets the x coordinate of the Lattice Point
	 *
	 * @return    The x value
	 */
	public int getX() {
		return this.x;
	}


	/**
	 *  Gets the y coordinate of the LatticePoint object
	 *
	 * @return    The y value
	 */
	public int getY() {
		return this.y;
	}


	/**
	 *  Gets the state of the Lattice Point
	 *
	 * @return    The state value
	 */
	public int getState() {
		return this.state;
	}


	/**
	 *  Check whether this point has changed
	 *
	 * @return    true if point has changed, false if not
	 */
	public boolean hasChanged() {
		return this.changed;
	}


	/**
	 *  Sets the state attribute of the LatticePoint object
	 *
	 * @param  state  The new state value
	 */
	public void setState(int state) {
		this.state = state;
	}


	/**
	 *  Set the boolean describing if the point has changed
	 *
	 * @param  changed  Whether the point has changed
	 */
	public void setChanged(boolean changed) {
		this.changed = changed;
	}
}

/**
 *  Class to hold a Cluster of Lattice Points
 *  Used for the Canvas to draw Wolff clusters
 *
 * @author    Richard Hanes
 */
class LatticeCluster {


	int oldState, newState;
	boolean[][] cluster;


	/**
	 *Constructor for the LatticeCluster object
	 *
	 * @param  oldState  The old state of the cluster
	 * @param  newState  The new state of the cluster
	 * @param  cluster   The boolean array holding the cluster information
	 */
	public LatticeCluster(int oldState, int newState, boolean[][] cluster) {
		this.oldState = oldState;
		this.newState = newState;
		this.cluster = cluster;
	}


	/**
	 *  Gets the oldState attribute of the LatticeCluster object
	 *
	 * @return    The oldState value
	 */
	public int getOldState() {
		return oldState;
	}


	/**
	 *  Gets the newState attribute of the LatticeCluster object
	 *
	 * @return    The newState value
	 */
	public int getNewState() {
		return newState;
	}


	/**
	 *  Is the current point in the cluster?
	 *
	 * @param  x  Coordinate of Lattice Point
	 * @param  y  Coordinate of Lattice Point
	 * @return    Whether or not the point is in the cluster
	 */
	public boolean isInCluster(int x, int y) {
		return cluster[x][y];
	}
}

/**
 *  Description of the Class
 *
 * @author    Richard Hanes
 */
class LatticeFrame extends Frame implements ActionListener {


	Lattice richie;
	LatticeCanvas pottsDisplay;
	LatticeControlPanel controls;


	/**
	 *Constructor for the LatticeFrame object
	 */
	public LatticeFrame() {
		setTitle("Emma's Lattice Simulator");
		this.controls = new LatticeControlPanel();

		int size = Integer.parseInt(controls.size.getText());
		int states = Integer.parseInt(controls.states.getText());
		double temp = Double.parseDouble(controls.temp.getText());
		double tempGradient = Double.parseDouble(controls.tempGradient.getText());

		this.richie = new Lattice(size, states, temp, tempGradient);
		this.pottsDisplay = new LatticeCanvas(richie);

		this.add(pottsDisplay, BorderLayout.CENTER);
		this.add(controls, BorderLayout.SOUTH);

		// add action listerners so the program 'listens' to the buttons
		controls.reset.addActionListener(this);
		controls.startPause.addActionListener(this);
		controls.algorithm.addActionListener(this);
		controls.update.addActionListener(this);

		// code to close the simulation if the close button is pushed
		addWindowListener(
			new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					System.exit(0);
				}
			});
	}


	/**
	 *  Listen for user input
	 *
	 * @param  e  The user event
	 */
	public void actionPerformed(ActionEvent e) {
		try {
			// If the update button is pushed
			if (e.getSource() == controls.update) {

				int size;

				int states;

				// Get the lattice parameters from the control panel
				states = Integer.parseInt(controls.states.getText());
				size = Integer.parseInt(controls.size.getText());

				// If the size has changed, update the lattice object
				if (richie.size != size) {
					richie.newSize = size;
					richie.sizeChanged = true;
					controls.size.setText(String.valueOf(size));
				}

				// If the number of states has changed update the lattice
				if (richie.states != states) {
					richie.newStates = states;
					richie.statesChanged = true;
					controls.states.setText(String.valueOf(states));
				}

				richie.temp = Double.parseDouble(controls.temp.getText());
				richie.tempGradient = Double.parseDouble(controls.tempGradient.getText());
				richie.calculateExponentials();
			}

			// The reset button doubles up as the update button when the simulation is paused
			// The update allows the user to watch each algorithm, step by step
			if (e.getSource() == controls.reset) {
				if (richie.running == true) {
					richie.resetMe = true;
				} else {
					if (richie.metropolis == true) {
						LatticePoint drawMe = richie.runMetropolis();
						if (drawMe.hasChanged()) {
							pottsDisplay.drawPoint(drawMe);
						}
					} else {
						LatticeCluster drawMe = richie.runWolff();
						pottsDisplay.drawCluster(drawMe);
					}
				}
			}

			if (e.getSource() == controls.startPause) {
				if (richie.running == false) {
					controls.startPause.setLabel("Pause");
					controls.reset.setLabel("Reset");
					richie.running = true;

				} else if (richie.running == true) {
					controls.startPause.setLabel("Continue");
					controls.reset.setLabel("Step");
					richie.running = false;
				}
			}
			if (e.getSource() == controls.algorithm) {
				if (richie.metropolis == true) {
					richie.metropolis = false;
					controls.temp.setText(String.valueOf(richie.temp));
					controls.algorithm.setLabel("Metropolis");
				} else {
					richie.metropolis = true;
					controls.temp.setText(String.valueOf(richie.temp));
					controls.algorithm.setLabel("Wolff");

				}
			}
		} catch (Exception ex) {
			System.out.println("Error in the Input");
			System.exit(1);
		}
	}


	/**
	 *  Run the Simulation
	 */
	public void runSimulation() {

		int i = 0;
		while (true) {

			// Check all of the 'Changed' flags every iteration and if any are true, update the lattice
			if (richie.running == true) {
				if (richie.sizeChanged == true) {
					richie.changeLatticeSize();
					richie.sizeChanged = false;
					pottsDisplay.repaint();
				}

				if (richie.statesChanged == true) {
					richie.changeNoOfStates();
					richie.statesChanged = false;
				}

				if (richie.resetMe == true) {
					richie.resetSimulation();
					pottsDisplay.repaint();
					richie.resetMe = false;
				}

				if (richie.tempGradient != 0) {
					if ((richie.temp + richie.tempGradient) > 0) {
						richie.temp += richie.tempGradient;
					} else {
						richie.temp = 0;
					}
					// Update the temperature on the control panel every 10000 iterations
					if (i == 10000) {
						controls.temp.setText(String.valueOf(richie.temp));
						i = 0;
					}
					i++;
				}

				if (richie.metropolis == true) {
					LatticePoint drawMe = richie.runMetropolis();
					if (drawMe.hasChanged()) {
						pottsDisplay.drawPoint(drawMe);
					}
				} else {
					LatticeCluster drawMe = richie.runWolff();
					pottsDisplay.drawCluster(drawMe);
				}
			}
		}
	}
}

/**
 *  The Canvas where the Lattice is Drawn
 *
 * @author    Richard Hanes
 */
class LatticeCanvas extends Canvas {


	Lattice richie;
	int size, states, w, h, x, y;
	Graphics safe;


	/**
	 * Constructor for the LatticeCanvas object
	 *
	 * @param  richie  The lattice to be drawn
	 */
	public LatticeCanvas(Lattice richie) {
		// set the background colour of the Canvas
		setBackground(Color.black);
		this.richie = richie;
	}


	/**
	 *  Allows for faster repainting
	 *
	 * @param  g  Graphics object of the Canvas
	 */
	public void update(Graphics g) {
		paint(g);
	}


	/**
	 *  Draw an individual point of the lattice
	 *
	 * @param  drawMe  The point to be drawn
	 */
	public void drawPoint(LatticePoint drawMe) {
		safe = getGraphics();
		Dimension d = getSize();
		int x;
		int y;
		int w;
		int h;
		int state;

		size = richie.getSize();
		states = richie.getStates();

		// determine size of each spin particle
		x = drawMe.getX();
		y = drawMe.getY();
		state = drawMe.getState();
		w = d.width / size;
		h = d.height / size;

		// paint the lattice
		x = w * x;
		y = h * y;
		safe.setColor(colour(state, states));
		safe.fillOval(x, y, w, h);

	}


	/**
	 *  Draw a cluster of Points
	 *
	 * @param  drawMe  The cluster of points to be drawn
	 */
	public void drawCluster(LatticeCluster drawMe) {
		size = richie.getSize();
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (drawMe.isInCluster(i, j)) {
					LatticePoint singleSpin = new LatticePoint(i, j, drawMe.getNewState(), true);
					drawPoint(singleSpin);
				}
			}
		}
	}


	/**
	 *  Draw the entire Lattice
	 *
	 * @param  g  Graphics object of the Canvas
	 */
	public void paint(Graphics g) {
		// get current size of Gui
		Dimension d = getSize();
		g.setColor(Color.black);

		g.fillRect(0, 0, d.width, d.height);

		size = richie.getSize();
		states = richie.getStates();

		// determine size of each spin particle
		w = d.width / size;
		h = d.height / size;

		// paint the lattice
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				x = w * i;
				y = h * j;
				g.setColor(colour(richie.latticeArray[i][j], states));
				g.fillOval(x, y, w, h);
			}
		}
	}


	/**
	 *  This method uses the rgb color model
	 *  to set the colour depending on ratio
	 *  of currentState to totalStates.
	 *
	 * @param  i       The state of the current point
	 * @param  states  The total number of states
	 * @return         A color object
	 */
	private Color colour(int i, int states) {

		float r;
		float g;
		float b;
		Color spinColor;

		// Vary the colours as the number of states increases
		switch (i) {
				case 7:
					spinColor = Color.yellow;
					break;
				case 8:
					spinColor = Color.blue;
					break;
				case 9:
					spinColor = Color.orange;
					break;
				case 10:
					spinColor = Color.pink;
					break;
				default:

					g = (float) i / (float) states;
					r = 1.0f - g;
					b = 0.0f;
					spinColor = new Color(r, g, b);
					break;
		}
		return spinColor;
	}
}

/**
 *  The control Panel at the bottom of the display
 *
 * @author    Richard Hanes
 */
class LatticeControlPanel extends Panel {


	Button startPause, algorithm, reset, update;
	TextField temp, tempGradient, states, size;
	Label tempLabel, statesLabel, sizeLabel, tempGradientLabel;


	/**
	 *Constructor for the LatticeControlPanel object
	 */
	public LatticeControlPanel() {
		setBackground(Color.gray);
		setLayout(new GridLayout(2, 6));

		//Initialise controls
		startPause = new Button("Pause");
		algorithm = new Button("Metropolis");
		reset = new Button("Reset");
		update = new Button("Update");

		tempLabel = new Label("Temperature", Label.RIGHT);
		tempGradientLabel = new Label("Temp Gradient", Label.RIGHT);
		statesLabel = new Label("States", Label.RIGHT);
		sizeLabel = new Label("Size", Label.RIGHT);

		temp = new TextField("0.99497", 4);
		tempGradient = new TextField("0.0", 4);
		states = new TextField("3", 4);
		size = new TextField("128", 4);

		//Add controls to panel
		add(startPause);
		add(reset);
		add(sizeLabel);
		add(size);
		add(statesLabel);
		add(states);
		add(algorithm);
		add(update);
		add(tempLabel);
		add(temp);
		add(tempGradientLabel);
		add(tempGradient);
	}
}

/**
 *  The Main class in the Program
 *
 * @author    Richard Hanes
 */
public class LatticeSimulator {


	/**
	 *  The main program for the LatticeSimulator class
	 *
	 * @param  args             The command line arguments
	 * @exception  IOException  Description of the Exception
	 */
	public static void main(String[] args) throws IOException {

		LatticeFrame potts = new LatticeFrame();

		// set size of Gui
		potts.setSize(650, 727);

		// make Gui visible
		potts.setVisible(true);
		potts.runSimulation();
	}
}



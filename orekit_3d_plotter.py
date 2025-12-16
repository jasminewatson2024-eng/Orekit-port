import os
import zipfile
import urllib.request
import sys
import math

# Initialize JVM and Orekit data (essential first steps)
import orekit
import java

# Import Orekit classes
from org.orekit.data import DataProvidersManager, DirectoryCrawler, DataContext
from org.orekit.time import AbsoluteDate, TimeScalesFactory, DateComponents, TimeComponents
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.utils import Constants, IERSConventions
from org.orekit.utils import PVCoordinates # For position and velocity vectors

# NEW IMPORTS FOR 3D PLOTTING
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # Required for 3D projection
import numpy as np # For numerical operations, e.g., creating Earth sphere

# --- 1. Initialize Orekit Data ---
def initialize_orekit_data():
    """
    Initializes the Orekit data context.
    This function will directly set the Orekit data path.
    """
    print("Initializing Orekit data...")
    orekit.initVM()
    print("Orekit JVM initialized successfully.")

    orekit_data_path = "C:\\OrekitProject\\orekit_data" # <<-- ENSURE THIS PATH IS CORRECT

    try:
        data_context = DataContext.getDefault()
        dm = data_context.getDataProvidersManager()
        data_provider = DirectoryCrawler(java.io.File(orekit_data_path))
        dm.addProvider(data_provider)
        print(f"Orekit data configured from: {orekit_data_path}")
    except Exception as e:
        print(f"Error initializing Orekit data provider from local path: {e}")
        print(f"Please ensure the folder '{orekit_data_path}' exists and contains the extracted Orekit data files.")
        sys.exit(1)

# --- 2. Create Satellite from TLE ---
def create_satellite_from_tle(tle_line1, tle_line2):
    """
    Creates a TLEPropagator and TLE object from two TLE lines.
    """
    tle = TLE(tle_line1, tle_line2)
    propagator = TLEPropagator.selectExtrapolator(tle)
    print(f"\nSatellite '{tle.getSatelliteNumber()}' created from TLE.")
    return propagator, tle

# --- 3. Propagate and Collect 3D Position Data ---
def propagate_and_collect_3d_positions(propagator, start_time, duration_seconds, step_seconds, utc):
    """
    Propagates the satellite and collects its X, Y, Z position coordinates
    in the GCRF frame at specified time steps.
    """
    print(f"\nPropagating for {duration_seconds / 3600:.2f} hours to collect 3D position data...")
    print(f"Starting from {start_time.toString()} (UTC) with {step_seconds} second steps.")

    gcrf_frame = FramesFactory.getGCRF() # Earth-Centered Inertial Frame
    end_time = start_time.shiftedBy(float(duration_seconds))

    x_coords = []
    y_coords = []
    z_coords = []

    current_date = start_time
    while current_date.compareTo(end_time) <= 0:
        try:
            propagated_state = propagator.propagate(current_date)
            position_gcrf = propagated_state.getPVCoordinates(gcrf_frame).getPosition()

            x_coords.append(position_gcrf.getX())
            y_coords.append(position_gcrf.getY())
            z_coords.append(position_gcrf.getZ())

            current_date = current_date.shiftedBy(float(step_seconds))
        except Exception as e:
            print(f"Error during propagation at {current_date}: {e}")
            break

    print(f"Collected {len(x_coords)} 3D position points.")
    return x_coords, y_coords, z_coords

# --- 4. Plot 3D Orbit ---
def plot_3d_orbit(x, y, z, output_filename="iss_3d_orbit.png"):
    """
    Creates a 3D plot of the satellite's orbit and saves it as a PNG image.
    This version focuses on making the orbit clearly visible.
    """
    print(f"\nGenerating 3D orbit plot: {output_filename}...")

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the orbit path - ADJUSTED LINEWIDTH
    ax.plot(x, y, z, color='red', linewidth=3, label='ISS Orbit') # Changed linewidth from 8 to 3

    # Add a marker for the start of the orbit
    if x and y and z:
        ax.scatter(x[0], y[0], z[0], color='green', s=100, label='Orbit Start', marker='o')

    # --- MODIFIED EARTH REPRESENTATION ---
    # Instead of a large sphere, plot a small black dot at the origin for Earth's center
    ax.scatter(0, 0, 0, color='black', s=200, label='Earth Center', marker='o') # Larger black dot for Earth center

    # --- CRUCIAL FIX: Adjusting Plot Limits to Fit the Orbit Tightly ---
    # Find the max absolute coordinate value from the orbit data only
    max_orbit_coord = max(np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z))) * 1.05 # Add a small buffer
    
    ax.set_xlim([-max_orbit_coord, max_orbit_coord])
    ax.set_ylim([-max_orbit_coord, max_orbit_coord])
    ax.set_zlim([-max_orbit_coord, max_orbit_coord])

    # Set labels and title
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('ISS 3D Orbit (GCRF)')
    ax.legend()

    # Set equal aspect ratio for a true 3D view
    ax.set_box_aspect([1,1,1]) # This ensures the plot box is cubic

    # Save the plot
    try:
        plt.savefig(output_filename, dpi=300) # Save with high resolution
        print(f"3D orbit plot saved successfully to {output_filename}")
        print(f"You can find this file in your C:\\OrekitProject folder.")
        print(f"Open '{output_filename}' to view the plot.")
        plt.show() # Display the interactive plot window
    except Exception as e:
        print(f"Error saving 3D plot: {e}")

# --- 5. Main Execution Block ---
if __name__ == "__main__":
    # Initialize Orekit data
    initialize_orekit_data()

    # Define time scales
    utc = TimeScalesFactory.getUTC()

    # Define ISS TLE
    iss_tle_line1 = "1 25544U 98067A   23204.57797746  .00008440  00000+0  15337-3 0  9997"
    iss_tle_line2 = "2 25544  51.6416 195.3854 0005709 130.4077 232.0963 15.49504547408798"

    # Create satellite propagator from TLE
    iss_propagator, iss_tle_obj = create_satellite_from_tle(iss_tle_line1, iss_tle_line2)

    # Define propagation time range for 3D plot (e.g., 1.5 orbits for ISS, approx 90 minutes per orbit)
    # ISS orbital period is approx 92 minutes. Let's do 2 hours to see a few loops.
    start_time_3d = AbsoluteDate(AbsoluteDate.J2000_EPOCH, 0.0, utc) # Start from J2000 epoch
    duration_for_3d_plot = float(2 * 3600) # 2 hours in seconds
    step_for_3d_plot = 60.0 # 1 minute steps

    # Propagate and collect 3D position data
    x_coords, y_coords, z_coords = propagate_and_collect_3d_positions(
        iss_propagator, start_time_3d, duration_for_3d_plot, step_for_3d_plot, utc
    )

    # Plot the 3D orbit if data was collected
    if x_coords:
        plot_3d_orbit(x_coords, y_coords, z_coords, output_filename="iss_3d_orbit.png")
    else:
        print("\nNo 3D orbit data to plot.")
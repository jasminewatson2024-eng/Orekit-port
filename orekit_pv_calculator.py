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
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.utils import Constants, IERSConventions
from org.orekit.utils import PVCoordinates # For position and velocity vectors

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

# --- 3. Calculate Position and Velocity ---
def calculate_pv(propagator, target_date, utc):
    """
    Calculates the satellite's position and velocity (PV) coordinates
    in Earth-Centered Inertial (GCRF) and Earth-Centered Earth-Fixed (ITRF) frames.
    """
    print(f"\nCalculating position and velocity for {target_date.toString()} (UTC)...")

    # Propagate to the target date
    propagated_state = propagator.propagate(target_date)

    # Get PVCoordinates in GCRF (Earth-Centered Inertial) frame
    # GCRF is a common quasi-inertial frame used for orbital mechanics
    gcrf_frame = FramesFactory.getGCRF()
    pv_gcrf = propagated_state.getPVCoordinates(gcrf_frame)

    # Get PVCoordinates in ITRF (Earth-Centered Earth-Fixed) frame
    # ITRF rotates with the Earth, useful for ground-relative positions
    itrf_frame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    pv_itrf = propagated_state.getPVCoordinates(itrf_frame)

    # Extract position and velocity vectors
    pos_gcrf = pv_gcrf.getPosition()
    vel_gcrf = pv_gcrf.getVelocity()

    pos_itrf = pv_itrf.getPosition()
    vel_itrf = pv_itrf.getVelocity()

    # Convert ITRF position to Geodetic (Lat, Lon, Alt) for human readability
    earth_body = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                  Constants.WGS84_EARTH_FLATTENING,
                                  itrf_frame) # Use ITRF as the body's frame
    geodetic_point = earth_body.transform(pos_itrf, itrf_frame, target_date)

    print("\n--- Position and Velocity Results ---")

    print("\nInertial Frame (GCRF):")
    print(f"  Position (X, Y, Z): ({pos_gcrf.getX():.3f} m, {pos_gcrf.getY():.3f} m, {pos_gcrf.getZ():.3f} m)")
    print(f"  Velocity (Vx, Vy, Vz): ({vel_gcrf.getX():.3f} m/s, {vel_gcrf.getY():.3f} m/s, {vel_gcrf.getZ():.3f} m/s)")
    print(f"  Speed: {pv_gcrf.getVelocity().getNorm():.3f} m/s")

    print("\nEarth-Fixed Frame (ITRF):")
    print(f"  Position (X, Y, Z): ({pos_itrf.getX():.3f} m, {pos_itrf.getY():.3f} m, {pos_itrf.getZ():.3f} m)")
    print(f"  Velocity (Vx, Vy, Vz): ({vel_itrf.getX():.3f} m/s, {vel_itrf.getY():.3f} m/s, {vel_itrf.getZ():.3f} m/s)")
    print(f"  Speed: {pv_itrf.getVelocity().getNorm():.3f} m/s")

    print("\nGeodetic Coordinates (from ITRF position):")
    print(f"  Latitude: {math.degrees(geodetic_point.getLatitude()):.4f}°")
    print(f"  Longitude: {math.degrees(geodetic_point.getLongitude()):.4f}°")
    print(f"  Altitude: {geodetic_point.getAltitude():.2f} m")

# --- 4. Main Execution Block ---
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

    # --- User Input for Target Date and Time ---
    print("\n--- Enter Target Date and Time for Calculation ---")
    year = int(input("Enter Year (e.g., 2023): "))
    month = int(input("Enter Month (1-12): "))
    day = int(input("Enter Day (1-31): "))
    hour = int(input("Enter Hour (0-23): "))
    minute = int(input("Enter Minute (0-59): "))
    second = float(input("Enter Second (0.0-59.999): "))

    # Create AbsoluteDate object from user input
    try:
        target_date_components = DateComponents(year, month, day)
        target_time_components = TimeComponents(hour, minute, second)
        target_absolute_date = AbsoluteDate(target_date_components, target_time_components, utc)
        print(f"\nTarget Date and Time set to: {target_absolute_date.toString()} (UTC)")
    except Exception as e:
        print(f"Error creating target date: {e}")
        print("Please ensure date and time inputs are valid.")
        sys.exit(1)

    # Calculate and display position and velocity
    calculate_pv(iss_propagator, target_absolute_date, utc)

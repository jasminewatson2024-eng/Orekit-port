# orekit_pass_predictor.py

# --- 1. Import necessary Orekit, Java, and Python modules ---
import os
import zipfile
import urllib.request
import sys # For exiting the script gracefully
import math # For standard mathematical functions like radians conversion

# It's crucial to import the Java Virtual Machine (JVM) and start it first.
# This must be done before importing any specific Orekit classes.
import orekit
import java # Import the java module for accessing Java classes

# Import Orekit classes
from org.orekit.data import DataProvidersManager, DirectoryCrawler, DataContext
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.orbits import KeplerianOrbit
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.propagation.events.handlers import EventHandler, StopOnIncreasing
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.utils import Constants, IERSConventions

# NEW IMPORT FOR PLOTTING
import folium

# --- 2. Initialize Orekit Data ---
def initialize_orekit_data():
    """
    Initializes the Orekit data context.
    This function will directly set the Orekit data path.
    """
    print("Initializing Orekit data...")
    orekit.initVM() # Initialize JVM as the very first step
    print("Orekit JVM initialized successfully.")

    # Define the absolute path to your locally downloaded Orekit data
    # THIS IS THE CRITICAL LINE TO ENSURE IT USES YOUR LOCAL DATA
    orekit_data_path = "C:\\OrekitProject\\orekit_data" # <<-- ENSURE THIS PATH IS CORRECT

    try:
        # Get the default data context
        data_context = DataContext.getDefault()

        # Get the data providers manager from the context
        dm = data_context.getDataProvidersManager()

        # Create a directory crawler for your data path
        data_provider = DirectoryCrawler(java.io.File(orekit_data_path)) # Use java.io.File

        # Add the directory crawler as a data provider
        dm.addProvider(data_provider) # Using addProvider, which is the standard Java method

        print(f"Orekit data configured from: {orekit_data_path}")

    except Exception as e:
        print(f"Error initializing Orekit data provider from local path: {e}")
        print(f"Please ensure the folder '{orekit_data_path}' exists and contains the extracted Orekit data files.")
        sys.exit(1) # Exit if data setup fails

# --- 3. Create Ground Station ---
def create_ground_station(station_name, latitude_deg, longitude_deg, altitude_m):
    """
    Creates an Orekit GroundStation object.
    """
    # Define Earth model (WGS84 ellipsoid)
    # Use IERSConventions.IERS_2010 for consistency with modern TLEs
    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, # Changed WGS88 to WGS84
                             Constants.WGS84_EARTH_FLATTENING, # Changed WGS88 to WGS84
                             FramesFactory.getITRF(IERSConventions.IERS_2010, True))

    # Define the geodetic point for the ground station
    station_point = GeodeticPoint(math.radians(latitude_deg),
                                  math.radians(longitude_deg),
                                  altitude_m)

    # Create the ground station
    ground_station = TopocentricFrame(earth, station_point, station_name)

    print(f"\nGround Station '{station_name}' created at:")
    print(f"  Latitude: {latitude_deg}°")
    print(f"  Longitude: {longitude_deg}°")
    print(f"  Altitude: {altitude_m} m")

    return ground_station, earth # Return earth as well for potential future use

# --- 4. Create Satellite from TLE ---
def create_satellite_from_tle(tle_line1, tle_line2):
    """
    Creates a TLEPropagator and TLE object from two TLE lines.
    """
    # Create TLE object
    tle = TLE(tle_line1, tle_line2)

    # Create TLEPropagator
    propagator = TLEPropagator.selectExtrapolator(tle)

    print(f"\nSatellite '{tle.getSatelliteNumber()}' created from TLE.")

    return propagator, tle

# --- 5. Custom Event Handler (Kept for completeness, but not used in this script) ---
# This class is no longer used in the ground track plotter, but kept to avoid errors if referenced elsewhere
class PassEventHandler(EventHandler):
    """
    Custom event handler for detecting and logging satellite passes.
    """
    def __init__(self, station_name):
        super().__init__()
        self.station_name = station_name
        self.pass_start_time = None
        self.pass_data = []

    def eventOccurred(self, s, detector, increasing):
        return EventHandler.Action.CONTINUE

    def getPasses(self):
        return self.pass_data

# --- 6. Satellite Ground Track Plotter Function ---
def calculate_satellite_ground_track(propagator, start_time, end_time, earth_body, step_seconds=60.0):
    """
    Calculates the satellite's ground track (latitude, longitude, altitude)
    at specified time steps within a given time range.
    """
    print(f"\nCalculating ground track for satellite {propagator.getTLE().getSatelliteNumber()}...")
    print(f"From {start_time} to {end_time} (UTC) with {step_seconds} second steps.")

    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

    ground_track_points = []
    current_date = start_time

    while current_date.compareTo(end_time) <= 0:
        try:
            propagated_state = propagator.propagate(current_date)
            position_itrf = propagated_state.getPVCoordinates(itrf).getPosition()
            geodetic_point = earth_body.transform(position_itrf, itrf, current_date)

            ground_track_points.append({
                "time": current_date,
                "latitude_deg": math.degrees(geodetic_point.getLatitude()),
                "longitude_deg": math.degrees(geodetic_point.getLongitude()),
                "altitude_m": geodetic_point.getAltitude()
            })
            current_date = current_date.shiftedBy(step_seconds)

        except Exception as e:
            print(f"Error during propagation at {current_date}: {e}")
            break

    print(f"\nGround track calculation complete. {len(ground_track_points)} points generated.")

    return ground_track_points

# --- 7. Helper Function to Split Ground Track for Plotting ---
def split_ground_track_for_plotting(ground_track_data, threshold=100):
    """
    Splits the ground track into segments to avoid drawing lines across the map
    when crossing the +/- 180 longitude meridian.
    """
    segments = []
    current_segment = []

    if not ground_track_data:
        return []

    current_segment.append((ground_track_data[0]['latitude_deg'], ground_track_data[0]['longitude_deg']))

    for i in range(1, len(ground_track_data)):
        prev_lon = ground_track_data[i-1]['longitude_deg']
        current_lon = ground_track_data[i]['longitude_deg']

        # Check for a large jump in longitude, indicating a meridian crossing
        # A jump greater than 'threshold' degrees (e.g., 100 degrees) suggests a wrap-around
        if abs(current_lon - prev_lon) > threshold:
            # If a jump is detected, end the current segment
            if current_segment:
                segments.append(current_segment)
            # Start a new segment with the current point
            current_segment = [(ground_track_data[i]['latitude_deg'], ground_track_data[i]['longitude_deg'])]
        else:
            # Otherwise, continue the current segment
            current_segment.append((ground_track_data[i]['latitude_deg'], ground_track_data[i]['longitude_deg']))

    # Add the last segment if it's not empty
    if current_segment:
        segments.append(current_segment)

    return segments

# --- 8. Plotting Function ---
def plot_ground_track(ground_track_data, output_filename="iss_ground_track.html"):
    """
    Plots the satellite ground track data on an interactive Folium map and saves it as an HTML file.
    """
    print(f"\nGenerating interactive map: {output_filename}...")

    m = folium.Map(location=[0, 0], zoom_start=2, tiles="OpenStreetMap")

    # Split the ground track data into segments to handle longitude wrap-around
    track_segments = split_ground_track_for_plotting(ground_track_data)

    # Add each segment as a separate polyline
    for segment in track_segments:
        folium.PolyLine(
            locations=segment,
            color='red',
            weight=2.5,
            opacity=1
        ).add_to(m)

    # Add markers for start and end points
    if ground_track_data:
        start_point = ground_track_data[0]
        end_point = ground_track_data[-1]

        folium.Marker(
            location=[start_point['latitude_deg'], start_point['longitude_deg']],
            popup=f"Start: {start_point['time'].toString()}",
            icon=folium.Icon(color='green', icon='play', prefix='fa')
        ).add_to(m)

        folium.Marker(
            location=[end_point['latitude_deg'], end_point['longitude_deg']],
            popup=f"End: {end_point['time'].toString()}",
            icon=folium.Icon(color='blue', icon='stop', prefix='fa')
        ).add_to(m)

    # Save the map to an HTML file
    try:
        m.save(output_filename)
        print(f"Map saved successfully to {output_filename}")
        print(f"You can find this file in your C:\\OrekitProject folder.")
        print(f"Open '{output_filename}' in any web browser to view the ground track.")
    except Exception as e:
        print(f"Error saving map: {e}")

# --- 9. Main Execution Block ---
if __name__ == "__main__":
    # Initialize Orekit data
    initialize_orekit_data()

    # Define time scales
    utc = TimeScalesFactory.getUTC()

    # Define Ground Station Parameters (not directly used in ground track plotting, but kept for context)
    austin_latitude_deg = 30.2849
    austin_longitude_deg = -97.7341
    austin_altitude_m = 150.00

    austin_station, earth_body = create_ground_station("UT Austin", austin_latitude_deg, austin_longitude_deg, austin_altitude_m)

    # Define ISS TLE
    iss_tle_line1 = "1 25544U 98067A   23204.57797746  .00008440  00000+0  15337-3 0  9997"
    iss_tle_line2 = "2 25544  51.6416 195.3854 0005709 130.4077 232.0963 15.49504547408798"

    # Create satellite propagator from TLE
    iss_propagator, iss_tle_obj = create_satellite_from_tle(iss_tle_line1, iss_tle_line2)

    # Define prediction time range (e.g., 2 days from now)
    current_time = AbsoluteDate(AbsoluteDate.J2000_EPOCH, 0.0, utc)
    duration_seconds = float(2 * 24 * 3600) # 2 days in seconds
    end_time = current_time.shiftedBy(duration_seconds)

    # Calculate the ground track data
    ground_track_data = calculate_satellite_ground_track(iss_propagator, current_time, end_time, earth_body, step_seconds=300.0)

    # Plot the ground track if data was generated
    if ground_track_data:
        plot_ground_track(ground_track_data, output_filename="iss_ground_track.html")
    else:
        print("\nNo ground track data to plot.")
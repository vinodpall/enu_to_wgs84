import numpy as np

# input
E, N, U = 0, 0, 0
lat_ref, lon_ref, h_ref = 31.6120280675, 120.785069774, 11.9676  # Latitude, longitude, and elevation of the reference point

a = 6378137.0
f = 1 / 298.257223563 
b = a * (1 - f)
e2 = 1 - (b**2 / a**2)

phi_ref = np.radians(lat_ref)
lambda_ref = np.radians(lon_ref)

N_ref = a / np.sqrt(1 - e2 * np.sin(phi_ref)**2)
X0 = (N_ref + h_ref) * np.cos(phi_ref) * np.cos(lambda_ref)
Y0 = (N_ref + h_ref) * np.cos(phi_ref) * np.sin(lambda_ref)
Z0 = (N_ref * (1 - e2) + h_ref) * np.sin(phi_ref)

R = np.array([
    [-np.sin(lambda_ref), np.cos(lambda_ref), 0],
    [-np.sin(phi_ref) * np.cos(lambda_ref), -np.sin(phi_ref) * np.sin(lambda_ref), np.cos(phi_ref)],
    [np.cos(phi_ref) * np.cos(lambda_ref), np.cos(phi_ref) * np.sin(lambda_ref), np.sin(phi_ref)]
])

dENU = np.array([E, N, U])
dECEF = np.dot(R.T, dENU)
X, Y, Z = X0 + dECEF[0], Y0 + dECEF[1], Z0 + dECEF[2]

lambda_target = np.arctan2(Y, X)
p = np.sqrt(X**2 + Y**2)
phi = np.arctan2(Z, p * (1 - e2))
N = a / np.sqrt(1 - e2 * np.sin(phi)**2)
h = p / np.cos(phi) - N

lat_target = np.degrees(phi)
lon_target = np.degrees(lambda_target)

print(f"WGS84：La = {lat_target:.8f}, Lo = {lon_target:.8f}, h = {h:.3f} m")

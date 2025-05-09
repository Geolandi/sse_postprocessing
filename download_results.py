import os
import requests
from requests.auth import HTTPBasicAuth
from bs4 import BeautifulSoup

# === Configuration ===
date = "2024-08-01"#2025-03-31"
BASE_URL = "https://near-real-time-sse.esc.cam.ac.uk/cascadia/" + date + "/"
USERNAME = "access"  # Replace with your username
PASSWORD = "zauSh~i3"  # Replace with your password
DOWNLOAD_FOLDER = "results/" + date
# make sure the folder exists
os.makedirs(DOWNLOAD_FOLDER, exist_ok=True)

# === STEP 1: Authenticate and Get File List ===
session = requests.Session()
response = session.get(BASE_URL, auth=HTTPBasicAuth(USERNAME, PASSWORD))

if response.status_code == 200:
    print("✅ Authentication successful!")
else:
    print(f"❌ Authentication failed! Status Code: {response.status_code}")
    print(response.text)
    exit()

# === STEP 2: Extract File Links ===
soup = BeautifulSoup(response.text, "html.parser")
file_links = [BASE_URL + a["href"] for a in soup.find_all("a", href=True) if not a["href"].startswith("?")]

# === STEP 3: Download Each File ===
os.makedirs(DOWNLOAD_FOLDER, exist_ok=True)

for file_url in file_links:
    filename = file_url.split("/")[-1]
    file_path = os.path.join(DOWNLOAD_FOLDER, filename)

    file_response = session.get(file_url, auth=HTTPBasicAuth(USERNAME, PASSWORD))
    if file_response.status_code == 200:
        with open(file_path, "wb") as f:
            f.write(file_response.content)
        print(f"✅ Downloaded {filename}")
    else:
        print(f"❌ Failed to download {filename} (Status: {file_response.status_code})")

print("🎉 All files downloaded successfully!")
